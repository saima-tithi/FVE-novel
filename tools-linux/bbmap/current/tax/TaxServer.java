package tax;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.net.InetSocketAddress;

import com.sun.net.httpserver.Headers;
import com.sun.net.httpserver.HttpExchange;
import com.sun.net.httpserver.HttpHandler;
import com.sun.net.httpserver.HttpServer;

import fileIO.ReadWrite;
import server.ServerTools;
import shared.Parser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import sketch.DisplayParams;
import sketch.Sketch;
import sketch.SketchMakerMini;
import sketch.SketchObject;
import sketch.SketchSearcher;
import sketch.SketchTool;
import sketch.Whitelist;
import stream.KillSwitch;
import structures.IntList;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.atomic.AtomicLong;

/**
 * @author Shijie Yao, Brian Bushnell
 * @date Dec 13, 2016
 *
 */
public class TaxServer {
	
	/*--------------------------------------------------------------*/
	/*----------------            Startup           ----------------*/
	/*--------------------------------------------------------------*/

	/** Command line entrance */
	public static void main(String[] args) throws Exception {
		Timer t=new Timer();
		TaxServer ts=new TaxServer(args);
		
		t.stop("Time: ");
		
		System.err.println("Ready!");
		
		//ts.begin();
	}
	
	/** Constructor */
	public TaxServer(String[] args) throws Exception {
		int port_=3068;
		String killCode_=null;
		
		TaxFilter.printNodesAdded=false;
		TaxFilter.REQUIRE_PRESENT=false; //Due to missing entries in TaxDump.
		
		//Create a parser object
		Parser parser=new Parser();
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b==null || b.equalsIgnoreCase("null")){b=null;}
			while(a.startsWith("-")){a=a.substring(1);} //Strip leading hyphens
			
			
			if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
			}else if(a.equals("verbose2")){
				verbose2=SketchObject.verbose2=Tools.parseBoolean(b);
			}else if(a.equals("table") || a.equals("gi") || a.equals("gitable")){
				tableFile=b;
				if("auto".equalsIgnoreCase(b)){tableFile=TaxTree.defaultTableFile();}
			}else if(a.equals("tree") || a.equals("taxtree")){
				treeFile=b;
				if("auto".equalsIgnoreCase(b)){treeFile=TaxTree.defaultTreeFile();}
			}else if(a.equals("accession")){
				accessionFile=b;
				if("auto".equalsIgnoreCase(b)){accessionFile=TaxTree.defaultAccessionFile();}
			}else if(a.equalsIgnoreCase("img")){
				imgFile=b;
			}else if(a.equals("reverse")){
				reverseOrder=Tools.parseBoolean(b);
			}else if(a.equals("domain")){
				domain=b;
				while(domain!=null && domain.endsWith("/")){domain=domain.substring(0, domain.length()-1);}
			}else if(a.equals("port")){
				port_=Integer.parseInt(b);
			}else if(a.equals("kill") || a.equals("killcode")){
				killCode_=b;
			}else if(a.equals("oldcode")){
				oldKillCode=b;
			}else if(a.equals("oldaddress")){
				oldAddress=b;
			}else if(a.equals("sketchonly")){
				sketchOnly=Tools.parseBoolean(b);
			}else if(a.equals("sketchthreads")){
				maxConcurrentSketchThreads=Integer.parseInt(b);
			}else if(a.equals("hashnames")){
				hashNames=Tools.parseBoolean(b);
			}else if(a.equals("printip")){
				printIP=Tools.parseBoolean(b);
			}else if(a.equals("countqueries")){
				countQueries=Tools.parseBoolean(b);
			}else if(a.equals("dbname")){
				SketchObject.defaultParams.dbName=b;
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else if(searcher.parse(arg, a, b, true)){
				//do nothing
			}else{
				throw new RuntimeException(arg);
			}
		}
		if("auto".equalsIgnoreCase(imgFile)){imgFile=TaxTree.defaultImgFile();}
		
		//Adjust SketchSearch rcomp and amino flags
		SketchObject.postParse();
		
		if(sketchOnly){
			hashNames=false;
			tableFile=null;
			accessionFile=null;
			imgFile=null;
		}
		
		port=port_;
		killCode=killCode_;
		
		//Fill some data objects
		USAGE=makeUsage();
		typeMap=makeTypeMap();
		commonMap=makeCommonMap();
		
		//Load the GI table
		if(tableFile!=null){
			outstream.println("Loading gi table.");
			GiToNcbi.initialize(tableFile);
		}
		
		//Load the taxTree
		if(treeFile!=null){
			tree=TaxTree.loadTaxTree(treeFile, outstream, hashNames);
			assert(tree.nameMap!=null || sketchOnly);
		}else{//The tree is required
			tree=null;
			throw new RuntimeException("No tree specified.");
		}
		//Set a default taxtree for sketch-related usage
		SketchObject.taxtree=tree;
		
		if(imgFile!=null){
			TaxTree.loadIMG(imgFile, false);
		}
		
		//Load accession files
		if(accessionFile!=null){
			AccessionToTaxid.tree=tree;
			outstream.println("Loading accession table.");
			AccessionToTaxid.load(accessionFile);
//			if(searcher.refFiles.isEmpty()){System.gc();}
		}
		
		//Load reference sketches
		hasSketches=searcher.refFileCount()>0;
		if(hasSketches){
			outstream.println("Loading sketches.");
			searcher.loadReferences(1, SketchObject.defaultParams.minEntropy);
//			System.gc();
		}
		
		{
			System.err.println("Clearing memory.");
			System.gc();
			Shared.printMemory();
		}
		
		//If there is a kill code, kill the old instance
		if(oldKillCode!=null && oldAddress!=null){
			killOldInstance();
		}
		
		//Wait for server initialization
		httpServer=initializeServer(2000, 7);
		assert(httpServer!=null);
		
		//Initialize handlers
		if(!sketchOnly){
			httpServer.createContext("/tax", new TaxHandler(false));
			httpServer.createContext("/stax", new TaxHandler(true));
			httpServer.createContext("/simpletax", new TaxHandler(true));
		}
		httpServer.createContext("/sketch", new SketchHandler());
		if(killCode!=null){
			httpServer.createContext("/kill", new KillHandler());
		}
//		httpServer.createContext("/help", new HelpHandler()); //Does not work as it is already caught by the default handler
		httpServer.createContext("/", new HelpHandler());
		httpServer.setExecutor(null); // creates a default executor
		
		//Start the server
		httpServer.start();
	}
	
	/** Kill a prior server instance */
	private void killOldInstance(){
		String result=null;
		try {
			result=ServerTools.sendAndReceive(oldKillCode.getBytes(), oldAddress);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.err.println("\nException suppressed; continuing.\n");
			return;
		}
		if(result==null || !result.equals("Success.")){
//			KillSwitch.kill("Bad kill result: "+result+"\nQuitting.\n");
			System.err.println("Bad kill result: "+result+"\nContinuing.\n");
		}
		ServerTools.pause(1000);
	}
	
	/** Iterative wait for server initialization */
	private HttpServer initializeServer(int millis0, int iterations){
		HttpServer server=null;
		InetSocketAddress isa=new InetSocketAddress(port);
		Exception ee=null;
		for(int i=0, millis=millis0; i<iterations && server==null; i++){
			try {
				server = HttpServer.create(isa, 0);
			} catch (java.net.BindException e) {//Expected
				System.err.println(e);
				System.err.println("\nWaiting "+millis+" ms");
				ee=e;
				ServerTools.pause(millis);
				millis=millis*2;
			} catch (IOException e) {//Not sure when this would occur...  it would be unexpected
				System.err.println(e);
				System.err.println("\nWaiting "+millis+" ms");
				ee=e;
				ServerTools.pause(millis);
				millis=millis*2;
			}
		}
		if(server==null){throw new RuntimeException(ee);}
		return server;
	}
	
	public void returnUsage(HttpExchange t){
		if(logUsage){System.err.println("usage");}
		ServerTools.reply(USAGE(), "text/plain", t, verbose2, 200);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Handlers           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Handles queries that fall through other handlers */
	class HelpHandler implements HttpHandler {
		
		@Override
		public void handle(HttpExchange t) throws IOException {
			returnUsage(t);
		}
		
	}
	
	/** Handles requests to kill the server */
	class KillHandler implements HttpHandler {
		
		@Override
		public void handle(HttpExchange t) throws IOException {
			
			//Parse the query from the URL
			String rparam=getRParam(t);
			InetSocketAddress remote=t.getRemoteAddress();
			
			if(testCode(t, rparam)){
				ServerTools.reply("Success.", "text/plain", t, verbose2, 200);
				System.err.println("Killed by remote address "+remote);
				//TODO: Perhaps try to close open resources such as the server
				KillSwitch.killSilent();
			}
			
			if(verbose){System.err.println("Bad kill from address "+remote);}
			ServerTools.reply(BAD_CODE, "text/plain", t, verbose2, 403);
		}
		
		/** Determines whether kill code was correct */
		private boolean testCode(HttpExchange t, String rparam){
			String[] params = rparam.split("/");
			if(verbose2){System.err.println(Arrays.toString(params));}
			
			if(killCode!=null){
				if(params.length>1){//URL mode
					return (params[1].equals(killCode));
				}else{//Body mode
					try {
						String code=ServerTools.receive(t);
						return (code!=null && code.equals(killCode));
					} catch (Exception e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
			return false;
		}
		
	}

	/** Listens for sketch comparison requests */
	class SketchHandler implements HttpHandler {
		
		@Override
		public void handle(HttpExchange t) throws IOException {
			
			final long startTime=System.nanoTime();
			
			//Parse the query from the URL
			String rparam=getRParam(t);
			
			assert(rparam.startsWith("sketch"));

			String body=getBody(t);
			if(!hasSketches){
				ServerTools.reply("\nERROR: This server has no sketches loaded.\n"
						+ "Please download the latest BBTools version to use SendSketch.\n", "text/plain", t, verbose2, 200);
				return;
			}
			
			if(rparam.startsWith("sketch/")){rparam=rparam.substring(7);}
			else{rparam=rparam.substring(6);}

			if(verbose2){
				System.err.println(rparam);
				System.err.println("rparam.startsWith(\"file/\"):"+rparam.startsWith("file/"));
			}
			
			
			//Toggle between local files and sketch transmission
			boolean fileMode=false;
			if(rparam.startsWith("file/")){
				if(verbose2){System.err.println("A");}
				rparam=rparam.substring(5);
				fileMode=true;
			}else{
				if(verbose2){System.err.println("B");}
			}
			incrementQueries(t, fileMode, false, false, -1);

			if(verbose2){System.err.println(rparam);}
			if(verbose2){System.err.println("fileMode="+fileMode);}
			
			//List of query sketches
			ArrayList<Sketch> sketches;
			DisplayParams params=SketchObject.defaultParams;
			
			if(verbose2){System.err.println("Found body: "+body);}
			if(body!=null && body.length()>0){
				if(fileMode && !body.startsWith("##")){
					body="##"+body;
				}
				try {
					params=SketchObject.defaultParams.parseDoubleHeader(body);
					if(verbose2){System.err.println("Passed parse params.");}
				} catch (Throwable e) {
					String s=Tools.toString(e);
					ServerTools.reply("\nERROR: \n"+ s,
							"text/plain", t, verbose2, 200);
					return;
				}
				if(!params.compatible()){
					ServerTools.reply("\nERROR: The sketch is not compatible with this server.\n"
							+ "Settings: k="+SketchObject.k+(SketchObject.k2>0 ? ","+SketchObject.k2 : "")+" amino="+SketchObject.amino+"\n"
							+ "You may need to download a newer version of BBTools; this is version "+Shared.BBMAP_VERSION_STRING,
							"text/plain", t, verbose2, 200);
					return;
				}
			}
			
			if(verbose2){System.err.println("Parsed params: "+params.toString());}
			
			if(fileMode){
				if(!new File(rparam).exists() && !rparam.startsWith("/")){
					String temp="/"+rparam;
					if(new File(temp).exists()){rparam=temp;}
				}
				sketches=loadSketchesFromFile(rparam, params);
			}else{
				sketches=loadSketchesFromBody(body);
			}
			if(verbose2){System.err.println("Loaded "+sketches.size()+" sketches.");}
			
			StringBuilder response=new StringBuilder();
			if(sketches==null || sketches.isEmpty()){
				response.append("Error.");
			}else{
				if(verbose2){
					System.err.println("Received "+sketches.get(0).name()+", size "+sketches.get(0).array.length);
					System.err.println("params: "+params);
				}
				searcher.compare(sketches, response, params, maxConcurrentSketchThreads); //This is where it gets stuck if comparing takes too long
				if(verbose2){System.err.println("Result: '"+response+"'");}
			}
			
			ServerTools.reply(response.toString(), "text/plain", t, verbose2, 200);
			t.close();
			
			final long stopTime=System.nanoTime();
			final long elapsed=stopTime-startTime;
			timeMeasurements.incrementAndGet();
			elapsedTime.addAndGet(elapsed);
			lastTime.set(elapsed);
		}
		
		private String getBody(HttpExchange t){
			InputStream is=t.getRequestBody();
			String s=ServerTools.readStream(is);
			return s;
		}
		
		private ArrayList<Sketch> loadSketchesFromBody(String body){
			//List of query sketches
			ArrayList<Sketch> sketches=null;
			
			if(body!=null && body.length()>0){
				sketches=searcher.loadSketchesFromString(body);
				if(Whitelist.exists()){
					for(Sketch sk : sketches){
						Whitelist.apply(sk);
					}
				}
			}
			return sketches;
		}
		
		private ArrayList<Sketch> loadSketchesFromFile(String fname, DisplayParams params){
			//List of query sketches
			ArrayList<Sketch> sketches=null;
			
			SketchTool tool=searcher.tool;
			if(tool.minKeyOccuranceCount!=params.minKeyOccuranceCount){tool=new SketchTool(SketchObject.targetSketchSize, params.minKeyOccuranceCount, params.printDepth);}
			if(verbose2){System.err.println("Loading sketches from file "+fname);}
			sketches=tool.loadSketches(fname, (SketchMakerMini)null, params.mode, params.samplerate, params.reads, params.minEntropy);
			
			return sketches;
		}
		
	}

	/** Handles taxonomy lookups */
	class TaxHandler implements HttpHandler {
		
		public TaxHandler(boolean skipNonCanonical_){
			skipNonCanonical=skipNonCanonical_;
		}
		
		@Override
		public void handle(HttpExchange t) throws IOException {
			final long startTime=System.nanoTime();
			
			//Parse the query from the URL
			String rparam=getRParam(t);

			String[] params = rparam.split("/");
			if(verbose2){System.err.println(Arrays.toString(params));}
			
			final String response=toResponse(skipNonCanonical, params, t);
			final String type=response.startsWith("{") ? "application/json" : "text/plain";
			
			ServerTools.reply(response, type, t, verbose2, 200);
			
			final long stopTime=System.nanoTime();
			final long elapsed=stopTime-startTime;
			timeMeasurements.incrementAndGet();
			elapsedTime.addAndGet(elapsed);
			lastTime.set(elapsed);
		}
		
		/** Only print nodes at canonical tax levels */
		public final boolean skipNonCanonical;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Helpers           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Parse the query from the URL */
	private String getRParam(HttpExchange t){
		String rparam = t.getRequestURI().toString();
		
		//Trim leading slashes
		while(rparam.startsWith("/")){
			rparam = rparam.substring(1);
		}
		
		//Trim trailing slashes
		while(rparam.endsWith("/")){
			rparam = rparam.substring(0, rparam.length()-1);
		}
		
		if(verbose){System.err.println(rparam);}
		return rparam;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------      Taxonomy Formatting     ----------------*/
	/*--------------------------------------------------------------*/
	
	/** All tax queries enter here from the handler */
	String toResponse(boolean skipNonCanonical, String[] params, HttpExchange t){
		if(params.length<3){
			if(params.length==2 && "advice".equalsIgnoreCase(params[1])){return TAX_ADVICE;}
			if(logUsage){System.err.println("usage");}
			return USAGE();
		}
		if(params.length>4){
			if(logUsage){System.err.println("usage");}
			return USAGE();
		}
		
		final String query=params[params.length-1];
		final String[] names=query.split(",");
		final boolean ancestor=(params.length>3 && params[2].equalsIgnoreCase("ancestor"));
//		System.err.println(params[2]+", "+ancestor);
		
		//Raw query type code
		final int type;
		{
			String typeS=params[1];
			Integer value=typeMap.get(typeS);
			if(value==null){
				if(typeS.equalsIgnoreCase("advice")){
					return TAX_ADVICE;
				}else{
					return "{\"error\": \"Bad type; should be gi, taxid, or name.\"}";
				}
			}
			type=value.intValue();
		}
		
		incrementQueries(t, false, skipNonCanonical, ancestor, type); //Ignores usage information.
		
		//Type code excluding formatting
		final int type2=type&15;
		final boolean plaintext=(type>=PT_OFFSET);
		final boolean semicolon=(type>=SC_OFFSET);
		
		if(verbose2){System.err.println("Type: "+type);}
		if(type2==NAME || type2==HEADER){
			for(int i=0; i<names.length; i++){
				names[i]=ServerTools.codeToSymbol(names[i]);
				if(type2==HEADER){
					if(names[i].startsWith("@") || names[i].startsWith(">")){names[i]=names[i].substring(1);}
				}
			}
			if(verbose2){System.err.println("Revised: "+Arrays.toString(names));}
		}
		
		if(ancestor){
			if(verbose2){System.err.println("toAncestor: "+Arrays.toString(names));}
			return toAncestor(type, names, plaintext, semicolon, query, skipNonCanonical, !skipNonCanonical);
		}
		
		if(semicolon){
			return toSemicolon(type, names, skipNonCanonical);
		}
		
		if(plaintext){
			return toText(type, names);
		}
		
		if(names.length==1){
			return toJson(type, names[0], skipNonCanonical, !skipNonCanonical).toString();
		}
		ArrayList<JsonObject> list=new ArrayList<JsonObject>();
		for(String name : names){
			list.add(toJson(type, name, skipNonCanonical, !skipNonCanonical));
		}
		return JsonObject.toString(list);
	}
	
	/** Look up common ancestor of terms */
	String toAncestor(final int type, final String[] names, boolean plaintext, boolean semicolon, String query, final boolean skipNonCanonical, boolean originalLevel){
		IntList ilist=toIntList(type, names);
		int id=FindAncestor.findAncestor(tree, ilist);
		TaxNode tn=(id>-1 ? tree.getNode(id) : null);
		if(tn==null){return new JsonObject(query, "error","Not found.").toString();}
		if(semicolon){
			return tree.toSemicolon(tn, skipNonCanonical);
		}
		if(plaintext){return ""+id;}
		
		JsonObject j=new JsonObject(query);
		j.add("name", tn.name);
		j.add("tax_id", ""+tn.id);
		j.add("level", ""+tn.levelStringExtended(originalLevel));
		while(tn!=null && tn.levelExtended!=TaxTree.LIFE_E && tn.id!=TaxTree.CELLULAR_ORGANISMS){
			if(!skipNonCanonical || tn.isSimple()){
				j.add(toJson(tn, originalLevel));
			}
			if(tn.pid==tn.id){break;}
			tn=tree.getNode(tn.pid);
		}
		return j.toString();
	}
	
	/** Format a reply as plaintext, comma-delimited, TaxID only */
	String toText(final int type, final String[] names){
		
		StringBuilder sb=new StringBuilder();
		String comma="";

		int type2=type&15;
		if(type2==GI){
			for(String name : names){
				sb.append(comma);
				TaxNode tn=getTaxNodeGi(Integer.parseInt(name));
				if(tn==null){sb.append("-1");}
				else{sb.append(tn.id);}
				comma=",";
			}
		}else if(type2==NAME){
			for(String name : names){
				sb.append(comma);
				TaxNode tn=getTaxNodeByName(name);
				if(tn==null){sb.append("-1");}
				else{sb.append(tn.id);}
				comma=",";
			}
		}else if(type2==NCBI){
			for(String name : names){
				sb.append(comma);
				TaxNode tn=getTaxNodeNcbi(Integer.parseInt(name));
				if(tn==null){sb.append("-1");}
				else{sb.append(tn.id);}
				comma=",";
			}
		}else if(type2==ACCESSION){
			for(String name : names){
				sb.append(comma);
				int ncbi=AccessionToTaxid.get(name);
				sb.append(ncbi);
				comma=",";
			}
		}else if(type2==HEADER){
			for(String name : names){
				sb.append(comma);
				TaxNode tn=getTaxNodeHeader(name);
				if(tn==null){sb.append("-1");}
				else{sb.append(tn.id);}
				comma=",";
			}
		}else if(type2==IMG){
			for(String name : names){
				sb.append(comma);
				int ncbi=TaxTree.imgToNcbi(Long.parseLong(name));
				sb.append(ncbi);
				comma=",";
			}
		}else{
			return "Bad type; should be pt_gi or pt_name; e.g. /tax/pt_gi/1234";
		}
		
		return sb.toString();
	}
	
	/** Format a reply as plaintext, semicolon-delimited, full lineage */
	String toSemicolon(final int type, final String[] names, boolean skipNonCanonical){
		
		StringBuilder sb=new StringBuilder();
		String comma="";
		
		int type2=type&15;
		if(type2==GI){
			for(String name : names){
				sb.append(comma);
				TaxNode tn=getTaxNodeGi(Integer.parseInt(name));
				if(tn==null){sb.append("Not found");}
				else{sb.append(tree.toSemicolon(tn, skipNonCanonical));}
				comma=",";
			}
		}else if(type2==NAME){
			for(String name : names){
				sb.append(comma);
				TaxNode tn=getTaxNodeByName(name);
				if(tn==null){sb.append("Not found");}
				else{sb.append(tree.toSemicolon(tn, skipNonCanonical));}
				comma=",";
			}
		}else if(type2==NCBI){
			for(String name : names){
				sb.append(comma);
				TaxNode tn=getTaxNodeNcbi(Integer.parseInt(name));
				if(tn==null){sb.append("Not found");}
				else{sb.append(tree.toSemicolon(tn, skipNonCanonical));}
				comma=",";
			}
		}else if(type2==ACCESSION){
			for(String name : names){
				sb.append(comma);
				final int tid=AccessionToTaxid.get(name);
				TaxNode tn=tree.getNode(tid, true);
				if(tn==null){sb.append("Not found");}
				else{sb.append(tree.toSemicolon(tn, skipNonCanonical));}
				comma=",";
			}
		}else if(type2==HEADER){
			for(String name : names){
				sb.append(comma);
				TaxNode tn=getTaxNodeHeader(name);
				if(tn==null){sb.append("Not found");}
				else{sb.append(tree.toSemicolon(tn, skipNonCanonical));}
				comma=",";
			}
		}else if(type2==IMG){
			for(String name : names){
				sb.append(comma);
				final int tid=TaxTree.imgToNcbi(Long.parseLong(name));
				TaxNode tn=tree.getNode(tid, true);
				if(tn==null){sb.append("Not found");}
				else{sb.append(tree.toSemicolon(tn, skipNonCanonical));}
				comma=",";
			}
		}else{
			return "Bad type; should be sc_gi or sc_name; e.g. /tax/sc_gi/1234";
		}
		
		return sb.toString();
	}
	
	/** Create a JsonObject from a String, including full lineage */
	JsonObject toJson(final int type, final String name, boolean skipNonCanonical, boolean originalLevel){
		TaxNode tn=null;
		
		if(type==GI){
			tn=getTaxNodeGi(Integer.parseInt(name));
		}else if(type==NAME){
			tn=getTaxNodeByName(name);
		}else if(type==NCBI){
			tn=getTaxNodeNcbi(Integer.parseInt(name));
		}else if(type==ACCESSION){
			int ncbi=AccessionToTaxid.get(name);
			tn=(ncbi>=0 ? tree.getNode(ncbi) : null);
		}else if(type==HEADER){
			tn=getTaxNodeHeader(name);
		}else if(type==IMG){
			final int tid=TaxTree.imgToNcbi(Long.parseLong(name));
			tn=tree.getNode(tid, true);
		}else{
			return new JsonObject(""+type,"error","Bad type; should be gi, taxid, or name; e.g. /tax/name/homo_sapiens");
		}
		if(verbose2){System.err.println("Got node: "+tn);}
		
		if(tn!=null){
			JsonObject j=new JsonObject(name);
			j.add("name", tn.name);
			j.add("tax_id", ""+tn.id);
			j.add("level", ""+tn.levelStringExtended(originalLevel));
			while(tn!=null && tn.levelExtended!=TaxTree.LIFE_E && tn.id!=TaxTree.CELLULAR_ORGANISMS){
//				System.err.println(tn+", "+(!skipNonCanonical)+", "+tn.isSimple());
				if(!skipNonCanonical || tn.isSimple()){
					j.add(toJson(tn, originalLevel));
//					System.err.println(j);
				}
				if(tn.pid==tn.id){break;}
				tn=tree.getNode(tn.pid);
			}
			return j;
		}
		return new JsonObject(name, "error","Not found.");
	}
	
	/** Create a JsonObject from a TaxNode, at that level only */
	JsonObject toJson(TaxNode tn, boolean originalLevel){
		JsonObject j=new JsonObject(tn.levelStringExtended(originalLevel));
		j.add("name",tn.name);
		j.add("tax_id",""+tn.id);
		return j;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Taxonomy Lookup       ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Convert a list of terms to a list of TaxIDs */
	IntList toIntList(final int type, final String[] names){
		IntList list=new IntList(names.length);
		int type2=type&15;
		if(type2==GI){
			for(String name : names){
				TaxNode tn=getTaxNodeGi(Integer.parseInt(name));
				if(tn!=null){list.add(tn.id);}
			}
		}else if(type2==NAME){
			for(String name : names){
				TaxNode tn=getTaxNodeByName(name);
				if(tn!=null){list.add(tn.id);}
			}
		}else if(type2==NCBI){
			for(String name : names){
				TaxNode tn=getTaxNodeNcbi(Integer.parseInt(name));
				if(tn!=null){list.add(tn.id);}
			}
		}else if(type2==ACCESSION){
			for(String name : names){
				int ncbi=AccessionToTaxid.get(name);
				if(ncbi>=0){list.add(ncbi);}
			}
		}else if(type2==IMG){
			for(String name : names){
				final int tid=TaxTree.imgToNcbi(Long.parseLong(name));
				if(tid>=0){list.add(tid);}
			}
		}else{
			throw new RuntimeException("{\"error\": \"Bad type\"}");
		}
		return list;
	}
	
	/** Look up a TaxNode by parsing the organism name */
	TaxNode getTaxNodeByName(String name){
		if(verbose2){System.err.println("Fetching node for "+name);}
		List<TaxNode> list=tree.getNodesByNameExtended(name);
		if(verbose2){System.err.println("Fetched "+list);}
		if(list==null){
			if(verbose2){System.err.println("Fetched in common map "+name);}
			String name2=commonMap.get(name);
			if(verbose2){System.err.println("Fetched "+name2);}
			if(name2!=null){list=tree.getNodesByName(name2);}
		}
		return list==null ? null : list.get(0);
	}
	
	/** Look up a TaxNode from the gi number */
	TaxNode getTaxNodeGi(int gi){
		int ncbi=-1;
		try {
			ncbi=GiToNcbi.getID(gi);
		} catch (Throwable e) {
			if(verbose){e.printStackTrace();}
		}
		return ncbi<0 ? null : getTaxNodeNcbi(ncbi);
	}
	
	//TODO:  Convert to using TaxTree.parseNodeFromHeader
	/** Look up a TaxNode by parsing the full header, assuming NCBI format */
	TaxNode getTaxNodeHeader(String header){
		return tree.parseNodeFromHeader(header, true);
	}
	
	/** Look up a TaxNode from the ncbi TaxID */
	TaxNode getTaxNodeNcbi(int ncbi){
		TaxNode tn=null;
		try {
			tn=tree.getNode(ncbi);
		} catch (Throwable e) {
			if(verbose){e.printStackTrace();}
		}
		return tn;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------     Data Initialization      ----------------*/
	/*--------------------------------------------------------------*/

	private static HashMap<String, Integer> makeTypeMap() {
		HashMap<String, Integer> map=new HashMap<String, Integer>(63);
		map.put("gi", GI);
		map.put("name", NAME);
		map.put("tax_id", NCBI);
		map.put("ncbi", NCBI);
		map.put("taxid", NCBI);
		map.put("id", NCBI);
		map.put("tid", NCBI);
		map.put("header", HEADER);
		map.put("accession", ACCESSION);
		map.put("img", IMG);
		map.put("pt_gi", PT_GI);
		map.put("pt_name", PT_NAME);
		map.put("pt_tax_id", PT_NCBI);
		map.put("pt_id", PT_NCBI);
		map.put("pt_tid", PT_NCBI);
		map.put("pt_ncbi", PT_NCBI);
		map.put("pt_taxid", PT_NCBI);
		map.put("pt_header", PT_HEADER);
		map.put("pt_header", PT_HEADER);
		map.put("pt_accession", PT_ACCESSION);
		map.put("pt_img", PT_IMG);
		map.put("sc_gi", SC_GI);
		map.put("sc_name", SC_NAME);
		map.put("sc_tax_id", SC_NCBI);
		map.put("sc_id", SC_NCBI);
		map.put("sc_tid", SC_NCBI);
		map.put("sc_ncbi", SC_NCBI);
		map.put("sc_taxid", SC_NCBI);
		map.put("sc_header", SC_HEADER);
		map.put("sc_header", SC_HEADER);
		map.put("sc_accession", SC_ACCESSION);
		map.put("sc_img", SC_IMG);
		
		return map;
	}
	
	public static HashMap<String, String> makeCommonMap(){
		HashMap<String, String> map=new HashMap<String, String>();
		map.put("human", "homo sapiens");
		map.put("cat", "felis catus");
		map.put("dog", "canis lupus familiaris");
		map.put("mouse", "mus musculus");
		map.put("cow", "bos taurus");
		map.put("bull", "bos taurus");
		map.put("horse", "Equus ferus");
		map.put("pig", "Sus scrofa domesticus");
		map.put("sheep", "Ovis aries");
		map.put("goat", "Capra aegagrus");
		map.put("turkey", "Meleagris gallopavo");
		map.put("fox", "Vulpes vulpes");
		map.put("chicken", "Gallus gallus domesticus");
		map.put("wolf", "canis lupus");
		map.put("fruitfly", "drosophila melanogaster");
		map.put("zebrafish", "Danio rerio");
		map.put("catfish", "Ictalurus punctatus");
		map.put("trout", "Oncorhynchus mykiss");
		map.put("salmon", "Salmo salar");
		map.put("tilapia", "Oreochromis niloticus");
		map.put("e coli", "Escherichia coli");
		map.put("e.coli", "Escherichia coli");

		map.put("lion", "Panthera leo");
		map.put("tiger", "Panthera tigris");
		map.put("bear", "Ursus arctos");
		map.put("deer", "Odocoileus virginianus");
		map.put("coyote", "Canis latrans");

		map.put("corn", "Zea mays subsp. mays");
		map.put("maize", "Zea mays subsp. mays");
		map.put("oat", "Avena sativa");
		map.put("wheat", "Triticum aestivum");
		map.put("rice", "Oryza sativa");
		map.put("potato", "Solanum tuberosum");
		map.put("barley", "Hordeum vulgare");
		map.put("poplar", "Populus alba");
		map.put("lettuce", "Lactuca sativa");
		map.put("beet", "Beta vulgaris");
		map.put("strawberry", "Fragaria x ananassa");
		map.put("orange", "Citrus sinensis");
		map.put("lemon", "Citrus limon");
		map.put("soy", "Glycine max");
		map.put("soybean", "Glycine max");
		map.put("grape", "Vitis vinifera");
		map.put("olive", "Olea europaea");
		map.put("cotton", "Gossypium hirsutum");
		map.put("apple", "Malus pumila");
		map.put("bannana", "Musa acuminata");
		map.put("tomato", "Solanum lycopersicum");
		map.put("sugarcane", "Saccharum officinarum");
		map.put("bean", "Phaseolus vulgaris");
		map.put("onion", "Allium cepa");
		map.put("garlic", "Allium sativum");
		
		map.put("pichu", "mus musculus");
		map.put("pikachu", "mus musculus");
		map.put("vulpix", "Vulpes vulpes");
		map.put("ninetails", "Vulpes vulpes");
		map.put("mareep", "Ovis aries");
		
		return map;
	}
	
	//Customize usage message to include domain
	private String makeUsage(){
		if(!sketchOnly){
			return "Welcome to the JGI taxonomy server!\n"
					+ "This service provides taxonomy information from NCBI taxID numbers, gi numbers, organism names, and accessions.\n"
					+ "The output is formatted as a Json object.\n"
					+ "Usage:\n\n"
					+ "All addresses are assumed to be prefixed by "+domain+", e.g.\n"
					+ domain+"/tax/name/homo_sapiens\n"
					+ "\n"
					+ "/tax/name/homo_sapiens will give taxonomy information for an organism name.\n"
					+ "Names are case-insensitive and underscores are equivalent to spaces.\n"
					+ "/tax/taxid/9606 will give taxonomy information for an NCBI taxID.\n"
					+ "/tax/gi/1234 will give taxonomy information from an NCBI gi number.\n"
					+ "/tax/accession/NZ_AAAA01000057.1 will give taxonomy information from an accession.\n"
					+ "/tax/header/ will accept an NCBI sequence header such as gi|7|emb|X51700.1| Bos taurus\n"
					+ "/tax/img/ will accept an IMG id such as 2724679250\n"
					+ "Vertical bars (|) cause problems for curl and can be replaced by tilde (~).\n"
					+ "\nComma-delimited lists are accepted for bulk queries, such as tax/gi/1234,7000,42\n"
					+ "For plaintext (non-Json) results, use the pt_ or sc_ prefix.\n"
					+ "pt_ will give just the taxID, while sc_ will give the whole lineage, semicolon-delimited.\n"
					+ "For example:\n\n"
					+ "/tax/pt_name/homo_sapiens\n"
					+ "/tax/sc_gi/1234\n"
					+ "\nTo find the common ancestor of multiple organisms, add /ancestor/. For example:\n"
					+ "/tax/taxid/ancestor/1234,5678,42\n"
					+ "/tax/name/ancestor/homo_sapiens,canis_lupus,bos_taurus\n"
					+ "\nFor a simplified taxonomic tree, use /simpletax or /stax instead of /tax.\n"
					+ "This will ignore unranked or uncommon levels like tribe and parvorder, and only display the following levels:\n"
					+ "SUBSPECIES, SPECIES, GENUS, FAMILY, ORDER, CLASS, PHYLUM, KINGDOM, SUPERKINGDOM, DOMAIN\n"
					+ "For example:\n"
					+ "/simpletax/taxid/ancestor/1234\n"
					+ "\nTo print taxonomy from the command line in Linux, use curl:\n"
					+ "curl http://taxonomy.jgi-psf.org/tax/taxid/9606\n"
					+ "\nLast restarted "+new Date()+"\n"
					+ "Running BBMap version "+Shared.BBMAP_VERSION_STRING+"\n";
		}else{
			StringBuilder sb=new StringBuilder();
			sb.append("Welcome to the JGI"+(SketchObject.defaultParams.dbName==null ? "" : " "+SketchObject.defaultParams.dbName)+" sketch server!\n");
//			if(dbName!=null){
//				sb.append("This server has the "+dbName+ " database loaded.\n");
//			}
			sb.append("\nUsage:\n\n");
			sb.append("sendsketch.sh in=file.fasta"+(SketchObject.defaultParams.dbName==null ? "" : " "+SketchObject.defaultParams.dbName.toLowerCase())+"\n\n");
			sb.append("SendSketch creates a sketch from a local sequence file, and sends the sketch to this server.\n");
			sb.append("The server receives the sketch, compares it to all sketches in memory, and returns the results.\n");
			sb.append("For files on the same system as the server, curl may be used with the absolute file path:\n\n");
			sb.append("curl "+domain+"/sketch/file/path\n");
			sb.append("e.g. curl "+domain+"/sketch/file//usr/data/assembly.fasta\n");
			sb.append("\nThis has the same effect as sendsketch.sh with the 'local' flag.\n");
			sb.append("\n");
			if(SketchObject.useWhitelist()){
				sb.append("This server is running in whitelist mode; for best results, use local queries.\n");
				sb.append("Remote queries should specify a larger-than-normal sketch size.\n\n");
			}else if(SketchObject.blacklist()!=null){
				sb.append("This server is running in blacklist mode, using "+new File(SketchObject.blacklist()).getName()+".\n\n");
			}
			sb.append("Last restarted "+new Date()+"\n");
			sb.append("Running BBMap version "+Shared.BBMAP_VERSION_STRING+"\n");
			sb.append("Settings: k="+SketchObject.k+(SketchObject.k2>0 ? ","+SketchObject.k2 : ""));
			if(SketchObject.amino){sb.append(" amino");}
			if(SketchObject.makeIndex){sb.append(" index");}
			if(SketchObject.useWhitelist()){sb.append(" whitelist");}
			if(SketchObject.blacklist()!=null){sb.append(" blacklist="+new File(SketchObject.blacklist()).getName());}
			sb.append('\n');
			return sb.toString();
		}
	}
	
	public String USAGE(){
		if(!countQueries){return USAGE;}
		final long uq=usageQueries.getAndIncrement();
		final long pt=plaintextQueries.get(), sc=semicolonQueries.get();
		final long i=internalQueries.get();
		final long l=localQueries.get();
		final long q=queries.get();
		final double avgTime=.000001*(elapsedTime.get()/(Tools.max(1.0, timeMeasurements.get())));//in milliseconds
		final double lastTimeD=.000001*lastTime.get();
		final long e=q-i;
		final long r=q-l;
		StringBuilder sb=new StringBuilder(USAGE.length()+500);
		sb.append(USAGE);

		sb.append('\n').append("Queries:   ").append(q);
		sb.append('\n').append("Usage:     ").append(uq);
		sb.append('\n').append("Avg time:  ").append(String.format("%.3f ms", avgTime));
		sb.append('\n').append("Last time: ").append(String.format("%.3f ms", lastTimeD));
		sb.append('\n');
		sb.append('\n').append("Internal:  ").append(i);
		sb.append('\n').append("External:  ").append(e);
		sb.append('\n');
		
		if(sketchOnly){
			sb.append('\n').append("Local:     ").append(l);
			sb.append('\n').append("Remote:    ").append(r);
		}else{
			sb.append('\n').append("gi:        ").append(giQueries.get());
			sb.append('\n').append("Name:      ").append(nameQueries.get());
			sb.append('\n').append("TaxID:     ").append(taxidQueries.get());
			sb.append('\n').append("Header:    ").append(headerQueries.get());
			sb.append('\n').append("Accession: ").append(accessionQueries.get());
			sb.append('\n').append("IMG:       ").append(imgQueries.get());
			sb.append('\n');
			sb.append('\n').append("Simple:    ").append(simpleQueries.get());
			sb.append('\n').append("Ancestor:  ").append(ancestorQueries.get());
			sb.append('\n');
			sb.append('\n').append("Json:      ").append(q-pt-sc);
			sb.append('\n').append("Plaintext: ").append(pt);
			sb.append('\n').append("Semicolon: ").append(sc);
		}
		sb.append('\n');
		return sb.toString();
	}
	
	public void incrementQueries(HttpExchange t, boolean local, boolean simple, boolean ancestor, int type){
		if(!countQueries){return;}
		queries.incrementAndGet();
		if(local){localQueries.incrementAndGet();}
		
		InetSocketAddress client=t.getRemoteAddress();
		InetSocketAddress server=t.getLocalAddress();

//		Headers reqh=t.getRequestHeaders();
//		Headers resh=t.getResponseHeaders();
//
//		System.err.println("\nRequest: ");
//		for(Entry<String, List<String>> entry : reqh.entrySet()){
//			System.err.println(entry.getKey()+" -> "+entry.getValue());
//		}
//		System.err.println("\nResponse: ");
//		for(Entry<String, List<String>> entry : resh.entrySet()){
//			System.err.println(entry.getKey()+" -> "+entry.getValue());
//		}
		
//		System.err.println(t.getRequestHeaders());
//		System.err.println(t.getResponseHeaders());
		
//		HttpContext context=t.getHttpContext();
//		context.getAttributes();
//		System.err.println(context.getAttributes());
		
//		t.getAttribute(arg0)
		
		if(printIP){System.err.println(client+"\t"+server);}
		
		if(type>=0){
			boolean plaintext=(type&PT_OFFSET)!=0;
			boolean semicolon=(type&SC_OFFSET)!=0;
			
			int type2=type&15;
			if(type2==GI){
				giQueries.incrementAndGet();
			}else if(type2==NAME){
				nameQueries.incrementAndGet();
			}else if(type2==NCBI){
				taxidQueries.incrementAndGet();
			}else if(type2==ACCESSION){
				accessionQueries.incrementAndGet();
			}else if(type2==IMG){
				imgQueries.incrementAndGet();
			}

			if(simple){simpleQueries.incrementAndGet();}
			if(ancestor){ancestorQueries.incrementAndGet();}
			
			if(plaintext){plaintextQueries.incrementAndGet();}
			if(semicolon){semicolonQueries.incrementAndGet();}
		}
		
		//This is for IPv4, class A.  Probably extends outside of Berkeley.
		String rs=client.toString();
		String ls=server.toString();
		
//		if(verbose2){
//			System.err.println(rs);
//			System.err.println(ls);
//		}
		
		boolean success=false;
		if(rs.contains("127.0.0.1")){
			Headers reqh=t.getRequestHeaders();
//			Headers resh=t.getResponseHeaders();
	
			String xff=reqh.getFirst("X-forwarded-for");
			if(xff!=null){
				if(xff.startsWith("128.")){
					internalQueries.incrementAndGet();
				}
				success=true;
			}
			
//			System.err.println("\nRequest: ");
//			for(Entry<String, List<String>> entry : reqh.entrySet()){
//				System.err.println(entry.getKey()+" -> "+entry.getValue());
//			}
		}
		
		if(!success){

			for(int i=0, max=Tools.max(rs.length(), ls.length()); i<max; i++){
				char rc=rs.charAt(i), lc=ls.charAt(i);
				if(rc!=lc){break;}
				if(rc=='.'){//IPv4
					internalQueries.incrementAndGet();
					break;
				}else if(rc==':'){//IPv6; probably depends on how long the mask is
					internalQueries.incrementAndGet();
					break;
				}
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------           Counters           ----------------*/
	/*--------------------------------------------------------------*/
	
	private AtomicLong queries=new AtomicLong(0);
	/** Same IP address mask */
	private AtomicLong internalQueries=new AtomicLong(0);
	/** Local filesystem sketch */
	private AtomicLong localQueries=new AtomicLong(0);

	private AtomicLong giQueries=new AtomicLong(0);
	private AtomicLong nameQueries=new AtomicLong(0);
	private AtomicLong taxidQueries=new AtomicLong(0);
	private AtomicLong headerQueries=new AtomicLong(0);
	private AtomicLong accessionQueries=new AtomicLong(0);
	private AtomicLong imgQueries=new AtomicLong(0);
	
	private AtomicLong plaintextQueries=new AtomicLong(0);
	private AtomicLong semicolonQueries=new AtomicLong(0);
	
	private AtomicLong simpleQueries=new AtomicLong(0);
	private AtomicLong ancestorQueries=new AtomicLong(0);

	private AtomicLong usageQueries=new AtomicLong(0);

	private AtomicLong elapsedTime=new AtomicLong(0);
	private AtomicLong timeMeasurements=new AtomicLong(0);
	private AtomicLong lastTime=new AtomicLong(0);
	
	/*--------------------------------------------------------------*/
	/*----------------            Params            ----------------*/
	/*--------------------------------------------------------------*/

	public boolean printIP=false;
	public boolean countQueries=true;

	/** Location of GiTable file */
	private String tableFile=null;
	/** Location of TaxTree file */
	private String treeFile=TaxTree.defaultTreeFile();
	/** Comma-delimited locations of Accession files */
	private String accessionFile=null;
	/** Location of IMG dump file */
	private String imgFile=null;
	
	/** Used for taxonomic tree traversal */
	private final TaxTree tree;
	
	/** Maps URL Strings to numeric query types */
	private final HashMap<String, Integer> typeMap;
	/** Maps common organism names to scientific names */
	private final HashMap<String, String> commonMap;
	
	/** Reverse order for tax lines */
	private boolean reverseOrder=true;
	
	/** Hash taxonomic names for lookup */
	private boolean hashNames=true;

	/** Kill code of prior server instance (optional) */
	private String oldKillCode=null;
	/** Address of prior server instance (optional) */
	private String oldAddress=null;

	/** Address of current server instance (optional) */
	public String domain="taxonomy.jgi-psf.org";
	
	public int maxConcurrentSketchThreads=16;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Listen on this port */
	public final int port;
	/** Code to validate kill requests */
	public final String killCode;
	
	public boolean sketchOnly=false;
	
	public final HttpServer httpServer;

	/** Bit to set for plaintext query types */
	public static final int PT_OFFSET=16;
	/** Bit to set for semicolon-delimited query types */
	public static final int SC_OFFSET=32;
	/** Request query types */
	public static final int UNKNOWN=0, GI=1, NAME=2, NCBI=3, HEADER=4, ACCESSION=5, IMG=6;
	/** Plaintext-response query types */
	public static final int PT_GI=GI+PT_OFFSET, PT_NAME=NAME+PT_OFFSET, PT_NCBI=NCBI+PT_OFFSET, 
			PT_HEADER=HEADER+PT_OFFSET, PT_ACCESSION=ACCESSION+PT_OFFSET, PT_IMG=IMG+PT_OFFSET;
	/** Semicolon-response query types */
	public static final int SC_GI=GI+SC_OFFSET, SC_NAME=NAME+SC_OFFSET, SC_NCBI=NCBI+SC_OFFSET, 
			SC_HEADER=HEADER+SC_OFFSET, SC_ACCESSION=ACCESSION+SC_OFFSET, SC_IMG=IMG+SC_OFFSET;
	
	/** Generic response when asking for tax advice */
	public static final String TAX_ADVICE="This site does not give tax advice.";
	/** Generic response for incorrect kill code */
	public static final String BAD_CODE="Incorrect code.";
	/** Generic response for badly-formatted queries */
	public final String USAGE;
	
	/** Tool for comparing query sketches to reference sketches */
	public final SketchSearcher searcher=new SketchSearcher();
	
	public final boolean hasSketches;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	private PrintStream outstream=System.err;
	/** Print verbose messages */
	public static boolean verbose=false, verbose2=false, logUsage=false;
	/** True if an error was encountered */
	public boolean errorState=false;
	
}
