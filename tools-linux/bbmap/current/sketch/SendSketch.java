package sketch;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import server.ServerTools;
import shared.Parser;
import shared.ReadStats;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ByteBuilder;

/**
 * Compares one or more input sketches to a set of reference sketches.
 * 
 * @author Brian Bushnell
 * @date July 29, 2016
 *
 */
public class SendSketch extends SketchObject {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		final boolean oldUnpigz=ReadWrite.USE_UNPIGZ;
		final int oldBufLen=Shared.bufferLen();
		
		//Create an instance of this class
		SendSketch ss=new SendSketch(args);
		
		//Run the object
		ss.process(t);
		
		ReadWrite.USE_UNPIGZ=oldUnpigz;
		Shared.setBufferLen(oldBufLen);
		assert(!ss.errorState) : "This program ended in an error state.";
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public SendSketch(String[] args){
		
		//Process any config files
		args=Parser.parseConfig(args);
		
		//Detect whether the uses needs help
		if(Parser.parseHelp(args, true)){
			printOptions();
			System.exit(0);
		}
		
		//Print the program name and arguments
//		outstream.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		
		//Set some shared static variables regarding PIGZ
		ReadWrite.USE_UNPIGZ=true;
		KILL_OK=true;
		
		//Create a parser object
		Parser parser=new Parser();
		parser.out1="stdout.txt";
		
		defaultParams.printFileName=true;
		boolean setBlacklist=false;
		boolean setLocal=false;
		boolean setPrintDepth=false;
		
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
			}else if(a.equals("in")){
				addFiles(b, in);
			}else if(a.equals("blacklist")){
				setBlacklist=true;
				parseSketchFlags(arg, a, b);
			}else if(parseSketchFlags(arg, a, b)){
				//Do nothing
			}else if(defaultParams.parse(arg, a, b)){
				//Do nothing
			}else if(a.equals("local")){
				local=Tools.parseBoolean(b);
				setLocal=true;
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Tools.parseKMG(b);
				//Set a variable here
			}else if(a.equals("address")){
				address=b;
				if(b.equalsIgnoreCase("nt")){
					address=ntAddress;
					blacklist=Blacklist.ntBlacklist();
				}else if(b.equalsIgnoreCase("refseq")){
					address=refseqAddress;
					blacklist=Blacklist.refseqBlacklist();
				}else if(b.equalsIgnoreCase("silva") || a.equalsIgnoreCase("ribo")){
					address=silvaAddress;
					blacklist=Blacklist.silvaBlacklist();
				}else if(b.equalsIgnoreCase("img")){
					address=imgAddress;
					blacklist=Blacklist.imgBlacklist();
				}
			}
			
			else if(a.equalsIgnoreCase("nt")){
				address=ntAddress;
				blacklist=Blacklist.ntBlacklist();
			}else if(a.equalsIgnoreCase("silva") || a.equalsIgnoreCase("ribo")){
				address=silvaAddress;
				blacklist=Blacklist.silvaBlacklist();
			}else if(a.equalsIgnoreCase("refseq")){
				address=refseqAddress;
				blacklist=Blacklist.refseqBlacklist();
			}else if(a.equalsIgnoreCase("img")){
				address=imgAddress;
				blacklist=Blacklist.imgBlacklist();
			}else if(a.equalsIgnoreCase("nr")){
				address=nrAddress;
				blacklist=Blacklist.nrBlacklist();
				amino=true;
			}
			
			else if(a.equals("name") || a.equals("taxname")){
				outTaxName=b;
			}else if(a.equals("name0")){
				outName0=b;
			}else if(a.equals("fname")){
				outFname=b;
			}else if(a.equals("taxid") || a.equals("tid")){
				outTaxID=Integer.parseInt(b);
			}else if(a.equals("spid")){
				outSpid=Integer.parseInt(b);
			}else if(a.equals("imgid")){
				outImgID=Integer.parseInt(b);
			}else if((a.startsWith("meta_") || a.startsWith("mt_")) && b!=null){
				if(outMeta==null){outMeta=new ArrayList<String>();}
				outMeta.add(a.substring(5)+":"+b);
			}
			
			else if(a.equals("outsketch") || a.equals("outs") || a.equals("sketchout") || a.equals("sketch")){
				outSketch=b;
			}
			
			else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}
			
			else if(b==null && in.isEmpty() && new File(arg).exists()){
				in.add(arg);
			}
			
			else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		
		while(address!=null && address.endsWith("/")){
			address=address.substring(0,  address.length()-1);
		}
		
		if(address.equals(nrAddress)){
			amino=true;
			if(!setK){k=9; k2=5;}
//			defaultParams.dbName="nr";
			assert(false) : "Need to set K.";
		}else if(address.equals(ntAddress)){
			if(!setK){k=defaultK; k2=24;}
//			defaultParams.dbName="nt";
			if(blacklist==null && !setBlacklist){blacklist=Blacklist.ntBlacklist();}
		}else if(address.equals(refseqAddress)){
			if(!setK){k=defaultK; k2=24;}
//			defaultParams.dbName="RefSeq";
			if(blacklist==null && !setBlacklist){blacklist=Blacklist.refseqBlacklist();}
		}else if(address.equals(silvaAddress)){
			if(!setK){k=defaultK; k2=0;}
//			defaultParams.dbName="Silva";
			if(blacklist==null && !setBlacklist){blacklist=Blacklist.silvaBlacklist();}
		}else if(address.equals(imgAddress)){
			if(!setK){k=defaultK; k2=24;}
//			defaultParams.dbName="IMG";
			if(blacklist==null && !setBlacklist){blacklist=Blacklist.imgBlacklist();}
		}
		
//		if(blacklist==null && !setBlacklist){
//			if(address==ntAddress){
//				blacklist=Blacklist.ntBlacklist;
//			}else if(address==refseqAddress){
//				blacklist=Blacklist.refseqBlacklist;
//			}else if(address==silvaAddress){
//				blacklist=Blacklist.silvaBlacklist;
//			}else if(address==imgAddress){
//				blacklist=Blacklist.imgBlacklist;
//			}
//		}
		
		if(local){blacklist=null;}
		
		postParse();
		
		{//Process parser fields
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;

			out=parser.out1;
		}
		
		//Ensure there is an input file
		if(in.isEmpty()){
			printOptions();
			throw new RuntimeException("Error - at least one input file is required.");
		}
		
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		ffout=FileFormat.testOutput(out, FileFormat.TEXT, null, false, overwrite, append, false);
		if(!ffout.stdio() && !defaultParams.setColors){defaultParams.printColors=false;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out, outSketch)){
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+
					out+", "+outSketch+"\n");
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in.toArray(new String[0]))){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
		
		tool=new SketchTool(targetSketchSize, defaultParams.minKeyOccuranceCount, defaultParams.trackCounts());
		
//		assert(false) : defaultParams.toString()+"\n"+k+", "+amino+", "+HASH_VERSION;
		if(verbose){
			if(local){System.err.println("Running in local mode.");}
			if(useWhitelist){System.err.println("Using a whitelist.");}
			if(blacklist!=null){System.err.println("Using a blacklist.");}
		}
		
		defaultParams.postParse(false);
		allowMultithreadedFastq=(in.size()==1 && Shared.threads()>2);
		if(!allowMultithreadedFastq){Shared.capBufferLen(40);}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public void process(Timer t){
		if(local){processLocal(t);}
		else{processRemote(t);}
	}
	
	private void processRemote(Timer t){
		Timer ttotal=new Timer();
		
		t.start();
		inSketches=tool.loadSketches_MT(defaultParams.mode, defaultParams.samplerate, defaultParams.reads, defaultParams.minEntropy, in);
		final int numLoaded=(inSketches.size());
		t.stop();
		outstream.println("Loaded "+numLoaded+" sketch"+(numLoaded==1 ? "" : "es")+" in \t"+t);
		t.start();
//		System.err.println(inSketches.get(0));
		
		if(ffout==null){return;}
		TextStreamWriter tsw=new TextStreamWriter(ffout);
		tsw.start();
		if(defaultParams.format==3){tsw.println(defaultParams.header());}
		
		ByteBuilder bb=new ByteBuilder();
		
		int cntr=0;
		for(Sketch sk : inSketches){
			
			if(sk.taxID<1 || sk.taxID>=minFakeID || outTaxID>0){sk.taxID=outTaxID;}
			if(outSpid>0){sk.spid=outSpid;}
			if(outImgID>0){sk.imgID=outImgID;}
			if(outTaxName!=null){sk.setTaxName(outTaxName);}
			if(outFname!=null){sk.setFname(outFname);}
			if(outName0!=null){sk.setName0(outName0);}
			sk.setMeta(outMeta);
			
			if(bb.length==0){bb.append(defaultParams.toString());}
			sk.toBytes(bb);
			cntr++;
			if(cntr>=100 || bb.length>500000){ //Don't allow too much data in a single transaction
				if(verbose){System.err.println("Sending:\n"+bb);}
				byte[] message=bb.toBytes();
				bb.clear();
				try {
//					System.err.println("Sending to "+address+"\n"+message+"\n"); //123
					String result=ServerTools.sendAndReceive(message, address);
					errorState|=checkForError(result);
					tsw.print(result);
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				cntr=0;
			}
		}
		
		if(bb.length>0){
			if(verbose){System.err.println("Sending:\n"+bb);}
			byte[] message=bb.toBytes();
			bb.clear();
			try {
				String result=ServerTools.sendAndReceive(message, address);
				errorState|=checkForError(result);
				tsw.println(result);
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		tsw.println();
		tsw.poison();
		
//		System.err.println("sending "+bb.toString());
		
		if(outSketch!=null){
			ByteStreamWriter bsw=new ByteStreamWriter(outSketch, overwrite, append, true, FileFormat.SKETCH);
			bsw.start();
			for(Sketch sk : inSketches){
				sk.toBytes(bb);
				bsw.print(bb);
				bb.clear();
			}
			bsw.poisonAndWait();
			errorState|=bsw.errorState;
		}
		
		tsw.waitForFinish();
		errorState|=tsw.errorState;
		
		t.stop();
//		outstream.println("\nRan "+(inSketches.size()*refSketches.size())+" comparisons in \t"+t);
		ttotal.stop();
		outstream.println("Total Time: \t"+ttotal);
	}
	
	private boolean checkForError(String s){
		if(s==null){return false;}
		return s.contains("java.io.IOException: Server returned HTTP response code:");
	}
	
	private void processLocal(Timer t){
		Timer ttotal=new Timer();
		
		t.start();
		
		if(ffout==null){return;}
		TextStreamWriter tsw=new TextStreamWriter(ffout);
		tsw.start();
		
		final String message=defaultParams.toString();
		for(String fname : in){
			String address2=address+"/file/"+new File(fname).getAbsolutePath();
			
			if(verbose){System.err.println("Sending:\n"+message);}
			try {
//				System.err.println("Sending to "+address2+"\n"+message+"\n"); //123
				String result=ServerTools.sendAndReceive(message.getBytes(), address2);
				tsw.print(result);
			} catch (Exception e) {
				// TODO Auto-generated catch block
				if(!suppressErrors){e.printStackTrace();}
			}
		}
		tsw.println();
		
//		System.err.println("sending "+bb.toString());
		
		tsw.poisonAndWait();
		errorState|=tsw.errorState;
		
		t.stop();
//		outstream.println("\nRan "+(inSketches.size()*refSketches.size())+" comparisons in \t"+t);
		ttotal.stop();
		outstream.println("Total Time: \t"+ttotal);
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private boolean addFiles(String a, Collection<String> list){
		int initial=list.size();
		if(a==null){return false;}
		File f=null;
		if(a.indexOf(',')>=0){f=new File(a);}
		if(f==null || f.exists()){
			list.add(a);
		}else{
			for(String s : a.split(",")){
				list.add(s);
			}
		}
		return list.size()>initial;
	}
	
	/** This is called if the program runs with no parameters */
	private void printOptions(){
		throw new RuntimeException("Please run the shellscript with no arguments for usage information."); //TODO
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private ArrayList<String> in=new ArrayList<String>();
	
	private String out="stdout.txt";
	private String outSketch=null;
	
	private final SketchTool tool;
	
	private ArrayList<Sketch> inSketches;

	private String address=refseqAddress;
	private boolean local=false;
	
	/*Override metadata */
	private String outTaxName=null;
	private String outFname=null;
	private String outName0=null;
	private int outTaxID=-1;
	private long outSpid=-1;
	private long outImgID=-1;
	private ArrayList<String> outMeta=null;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary output file */
	private final FileFormat ffout;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	private PrintStream outstream=System.err;
	/** Print verbose messages */
	public static boolean verbose=false;
	/** True if an error was encountered */
	public boolean errorState=false;
	/** Overwrite existing output files */
	private boolean overwrite=false;
	/** Append to existing output files */
	private boolean append=false;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Don't print caught exceptions */
	public static boolean suppressErrors=false;
	
	private static final String ntAddress="https://nt-sketch.jgi-psf.org/sketch"; //port 3071
	private static final String refseqAddress="https://refseq-sketch.jgi-psf.org/sketch"; //port 3072
	private static final String silvaAddress="https://ribo-sketch.jgi-psf.org/sketch"; //port 3073
	private static final String imgAddress="https://img-sketch.jgi-psf.org/sketch";
	private static final String nrAddress="https://nr-sketch.jgi-psf.org/sketch";
	
}
