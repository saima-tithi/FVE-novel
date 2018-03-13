package sketch;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedHashSet;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;

import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import shared.Parser;
import shared.ReadStats;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import tax.TaxTree;

/**
 * Compares one or more input sketches to a set of reference sketches.
 * 
 * @author Brian Bushnell
 * @date July 29, 2016
 *
 */
public class CompareSketch extends SketchObject {
	

	
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
		
		final int oldBufLen=Shared.bufferLen();
		
		//Create an instance of this class
		CompareSketch cs=new CompareSketch(args);
		
		//Run the object
		cs.process(t);

		Shared.setBufferLen(oldBufLen);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public CompareSketch(String[] args){
		
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
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=Shared.threads();
		KILL_OK=true;
		
		//Create a parser object
		Parser parser=new Parser();
		parser.out1="stdout.txt";
		
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
			}else if(parseSketchFlags(arg, a, b)){
				//Do nothing
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Tools.parseKMG(b);
				//Set a variable here
			}else if(a.equals("ordered")){
				ordered=Tools.parseBoolean(b);
			}else if(a.equals("alltoall") || a.equals("ata")){
				allToAll=Tools.parseBoolean(b);
			}else if(a.equals("compareself") || a.equals("includeself")){
				compareSelf=Tools.parseBoolean(b);
			}else if(a.equals("printmemory")){
				printMemory=Tools.parseBoolean(b);
			}
			
			else if(a.equals("taxtree") || a.equals("tree")){
				taxTreeFile=b;
				if("auto".equalsIgnoreCase(b)){taxTreeFile=TaxTree.defaultTreeFile();}
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
			
			else if(searcher.parse(arg, a, b, false)){
				//do nothing
			}
			
			else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}
			
			else if(searcher.parse(arg, a, b, true)){
				//do nothing
			}
			
			else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		SketchObject.postParse();
		
		if(in.isEmpty() && args.length>0 && !allToAll){ //Allows first argument to be used as the input file without in= flag
			String x=args[0];
			if(x.indexOf('=')<0 && new File(x).exists() && searcher.refFiles.contains(x)){
				searcher.refFiles.remove(x);
				in.add(x);
			}
		}
		
		{//Process parser fields
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;

			out=parser.out1;
		}
		
		if(allToAll){
			 LinkedHashSet<String> set=new LinkedHashSet<String>();
			 set.addAll(in);
			 set.addAll(searcher.refFiles);
			 in.clear();
			 searcher.refFiles.clear();
			 in.addAll(set);
			 searcher.refFiles.addAll(set);
		}
		
		//Ensure there is an input file
		if(in.isEmpty()){
			printOptions();
			throw new RuntimeException("Error - at least one input file is required.");
		}
		
		//Ensure there is an ref file
		if(searcher.refFiles.isEmpty()){
			printOptions();
			throw new RuntimeException("Error - at least one reference file is required.");
		}
		
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		ffout=FileFormat.testOutput(out, FileFormat.TEXT, null, false, overwrite, append, false);
		if(!ffout.stdio() && !defaultParams.setColors){defaultParams.printColors=false;}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, taxTreeFile)){
			throw new RuntimeException("\nCan't read to some input files.\n");
		}
		if(!Tools.testInputFiles(true, false, in.toArray(new String[0]))){
			if(in.size()==1){
				String s=in.get(0);
				String s1=s.replaceFirst("#", "1"), s2=s.replaceFirst("#", "2");
				Tools.testInputFiles(true, false, s1, s2);
			}else{
				throw new RuntimeException("\nCan't read to some input files.\n");
			}
		}
		
//		assert(makeIndex || defaultParams.printContam2==false) : "Contam2 requires the flag index=t";
		
		if(taxTreeFile!=null){setTaxtree(taxTreeFile);}
		defaultParams.postParse(true);
		allowMultithreadedFastq=in.size()<2 && !allToAll;
		if(!allowMultithreadedFastq){Shared.capBufferLen(40);}
		
//		assert(defaultParams.checkValid());
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public void process(Timer t){
		Timer ttotal=new Timer();
		
		t.start();
		
		outstream.println("Loading sketches.");
		searcher.makeTool(1, false);
		SketchTool tool=new SketchTool(targetSketchSize, defaultParams.minKeyOccuranceCount, defaultParams.trackCounts());
		
		if(!useWhitelist || allToAll){
			if(allToAll){
				makeIndex=makeIndex || defaultParams.needIndex() || searcher.autoIndex;
				searcher.loadReferences(defaultParams.minKeyOccuranceCount, defaultParams.minEntropy);
				inSketches=(ArrayList<Sketch>) searcher.refSketches.clone();
			}else{
				inSketches=tool.loadSketches_MT(defaultParams.mode, defaultParams.samplerate, defaultParams.reads, defaultParams.minEntropy, in);
				
				for(Sketch sk : inSketches){
					if(sk.taxID<1 || sk.taxID>=minFakeID || outTaxID>0){sk.taxID=outTaxID;}
					if(outSpid>0){sk.spid=outSpid;}
					if(outImgID>0){sk.imgID=outImgID;}
					if(outTaxName!=null){sk.setTaxName(outTaxName);}
					if(outFname!=null){sk.setFname(outFname);}
					if(outName0!=null){sk.setName0(outName0);}
					sk.setMeta(outMeta);
				}
				
				if(outTaxID>0){
					for(Sketch sk : inSketches){
						if(sk.taxID<1 || sk.taxID>=minFakeID){sk.taxID=outTaxID;}
					}
				}
				makeIndex=(searcher.autoIndex && inSketches.size()>8) || defaultParams.needIndex();
				searcher.loadReferences(defaultParams.minKeyOccuranceCount, defaultParams.minEntropy);
			}
		}else{
			//assert(searcher.makeIndex && !searcher.autoIndex) : "whitelist=t requires index=t";
			makeIndex=true;
			searcher.loadReferences(defaultParams.minKeyOccuranceCount, defaultParams.minEntropy);
			inSketches=tool.loadSketches_MT(defaultParams.mode, defaultParams.samplerate, defaultParams.reads, defaultParams.minEntropy, in);
		}
		
		final int numLoaded=(inSketches.size()+searcher.refSketches.size())/(allToAll ? 2 : 1);
		t.stop();
		outstream.println("Loaded "+numLoaded+" sketch"+(numLoaded==1 ? "" : "es")+" in \t"+t);
		if(printMemory){
			System.gc();
			Shared.printMemory();
		}
		t.start();

		
		TextStreamWriter tsw=(ffout==null ? null : new TextStreamWriter(ffout));
		if(tsw!=null){tsw.start();}
		if(defaultParams.format==3){tsw.println(defaultParams.header());}

		boolean success=true;
		final int inSize=inSketches.size();
		if(inSize==1 || Shared.threads()<2 || inSize<Shared.threads()){
			StringBuilder sb=new StringBuilder();
			success=searcher.compare(inSketches, sb, defaultParams, Shared.threads());
			success&=(!searcher.errorState);
			if(tsw!=null){tsw.println(sb);}
		}else{//More sketches than threads, and more than one thread
			final int threads=Tools.min(Shared.threads(), inSize);
			
			ArrayList<CompareThread> alct=new ArrayList<CompareThread>(threads);
			AtomicInteger next=new AtomicInteger(0);
			for(int i=0; i<threads; i++){
				alct.add(new CompareThread(i, next, tsw));
			}
			for(CompareThread ct : alct){ct.start();}
			for(CompareThread ct : alct){

				//Wait until this thread has terminated
				while(ct.getState()!=Thread.State.TERMINATED){
					try {
						//Attempt a join operation
						ct.join();
					} catch (InterruptedException e) {
						e.printStackTrace();
					}
				}

				synchronized(ct){
					success&=ct.success;
				}
			}
			alct=null;
		}
		
		//Track whether any threads failed
		if(!success){errorState=true;}
		if(tsw!=null){tsw.poisonAndWait();}
		errorState|=tsw.errorState;
		
		t.stop();
		long comparisons=(makeIndex ? searcher.comparisons.get() : inSketches.size()*(long)searcher.refSketches.size());
		outstream.println("\nRan "+comparisons+" comparison"+(comparisons==1 ? "" : "s")+" in \t"+t);
		ttotal.stop();
		outstream.println("Total Time: \t"+ttotal);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private void addFiles(String a, Collection<String> list){
		if(a==null){return;}
		File f=null;
		if(a.indexOf(',')>=0){f=new File(a);}
		if(f==null || f.exists()){
			list.add(a);
		}else{
			for(String s : a.split(",")){list.add(s);}
		}
	}
	
	/** This is called if the program runs with no parameters */
	private void printOptions(){
		throw new RuntimeException("Please run the shellscript with no arguments for usage information."); //TODO
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	private class CompareThread extends Thread {
		
		CompareThread(final int tid_, final AtomicInteger nextSketch_, TextStreamWriter tsw_){
			tid=tid_;
			nextSketch=nextSketch_;
			tsw=tsw_;
		}
		
		public void run(){
			success=false;
			final int inLim=inSketches.size();
			
			for(int inNum=nextSketch.getAndIncrement(); inNum<inLim; inNum=nextSketch.getAndIncrement()){
				Sketch a=inSketches.get(inNum);
				assert(buffer.cbs==null); //Because this sketch will only be used by one thread at a time, so per-buffer bitsets are not needed. 
				SketchResults sr=searcher.processSketch(a, buffer, fakeID, map, defaultParams, 1);
				
				if(tsw!=null){
					synchronized(tsw){
						if(ordered){
							StringBuilder sb=sr.toText(defaultParams);
							tsw.writeOrdered(sb, inNum);
						}else{
							sr.writeResults(defaultParams, tsw);
						}
					}
				}
			}
			synchronized(this){success=true;}
		}
		
		private final int tid;
		private final CompareBuffer buffer=new CompareBuffer(false);

		private final AtomicInteger nextSketch;
		private final AtomicInteger fakeID=new AtomicInteger(minFakeID);
		private ConcurrentHashMap<Integer, Comparison> map=new ConcurrentHashMap<Integer, Comparison>(101);
		final TextStreamWriter tsw;
		
		boolean success=false;
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private ArrayList<String> in=new ArrayList<String>();
	
	private String out="stdout.txt";
	
	private String taxTreeFile=null;
	
	private ArrayList<Sketch> inSketches;
	
	public final SketchSearcher searcher=new SketchSearcher();
	
	private boolean printMemory=false;
	
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
	private boolean ordered=true;
	
}
