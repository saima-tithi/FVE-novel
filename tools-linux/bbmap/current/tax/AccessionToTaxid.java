package tax;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicIntegerArray;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import kmer.HashBuffer;
import kmer.KmerTableSet;
import shared.Parser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentGenericReadInputStream;
import stream.FastaReadInputStream;

/**
 * @author Brian Bushnell
 * @date December 16, 2016
 *
 */
public class AccessionToTaxid {
	
	public static void load(String files){
		main(new String[] {"in="+files});
	}
	
	public static void main(String[] args){
		Timer t=new Timer();
		AccessionToTaxid sample=new AccessionToTaxid(args);
		sample.process(t);
	}
	
	public AccessionToTaxid(String[] args){
		
		args=Parser.parseConfig(args);
		if(Parser.parseHelp(args, true)){
			printOptions();
			System.exit(0);
		}
		
		outstream.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=Shared.threads();
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b==null || b.equalsIgnoreCase("null")){b=null;}
			while(a.startsWith("-")){a=a.substring(1);} //In case people use hyphens

			if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("stripunderscore")){
//				STRIP_UNDERSCORE=Tools.parseBoolean(b);
				assert(false) : "stripunderscore is disabled.";
			}else if(a.equals("usetables")){
//				USE_TABLES=Tools.parseBoolean(b);
			}else if(a.equals("countlengths")){
				boolean x=Tools.parseBoolean(b);
				lengthCounts=x ? new AtomicIntegerArray(20) : null;
			}else if(a.equals("in")){
				String[] temp=b.split(",");
				for(String s : temp){in.add(s);}
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else if(b==null){
				if(new File(arg).exists()){
					in.add(arg);
				}
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			overwrite=parser.overwrite;

			out=parser.out1;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in==null || in.size()==0){
			printOptions();
			throw new RuntimeException("Error - at least one input file is required.");
		}
		
		if(!ByteFile.FORCE_MODE_BF2){
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF1=true;
		}

		if(out!=null && out.equalsIgnoreCase("null")){out=null;}
		
		if(!Tools.testOutputFiles(overwrite, false, false, out)){
			outstream.println((out==null)+", "+out);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out+"\n");
		}

		ffout=FileFormat.testOutput(out, FileFormat.TXT, null, true, overwrite, false, false);
		ffin=new FileFormat[in.size()];
		for(int i=0; i<in.size(); i++){
			String s=in.get(i);
			if(!new File(s).exists()){
				if(s.startsWith("shrunk.") && new File(s.substring(7)).exists()){
					s=s.substring(7);
				}
			}
			ffin[i]=FileFormat.testInput(s, FileFormat.TXT, null, true, true);
		}
	}
	
	@SuppressWarnings("unchecked")
	void process(Timer t){

//		if(USE_MAPS){
			assert(maps==null);
			maps=new HashMap[128];
			for(int i=0; i<maps.length; i++){
				maps[i]=new HashMap<String, Integer>();
			}
//		}

		assert(tables==null);
		if(USE_TABLES){
			tables=new KmerTableSet(new String[] {"ways=31"}, 12);
			tables.allocateTables();
		}
		
		spawnThreads();
		
		t.stop();
		
		double rpnano=linesProcessed/(double)(t.elapsed);
		double bpnano=bytesProcessed/(double)(t.elapsed);

		String rpstring=(linesProcessed<100000 ? ""+linesProcessed : linesProcessed<100000000 ? (linesProcessed/1000)+"k" : (linesProcessed/1000000)+"m");
		String bpstring=(bytesProcessed<100000 ? ""+bytesProcessed : bytesProcessed<100000000 ? (bytesProcessed/1000)+"k" : (bytesProcessed/1000000)+"m");

		while(rpstring.length()<8){rpstring=" "+rpstring;}
		while(bpstring.length()<8){bpstring=" "+bpstring;}
		
		outstream.println("Time:                         \t"+t);
		outstream.println("Lines Processed:    "+rpstring+" \t"+String.format(Locale.ROOT, "%.2fk lines/sec", rpnano*1000000));
		outstream.println("Bytes Processed:    "+bpstring+" \t"+String.format(Locale.ROOT, "%.2fm bytes/sec", bpnano*1000));
		
		outstream.println();
		outstream.println("Valid Lines:       \t"+linesValid);
		outstream.println("Invalid Lines:     \t"+(linesProcessed-linesValid));

		if(lengthCounts!=null){
			outstream.println();
			outstream.println("Length counts:");

			for(int i=0; i<lengthCounts.length(); i++){
				int count=lengthCounts.get(i);
				if(count>0){outstream.println(i+"\t"+count);}
			}
		}

		if(symbolCounts!=null){
			outstream.println();
			outstream.println("Symbols:");
			
			String comma="";
			for(int i=0; i<symbolCounts.length(); i++){
				int count=symbolCounts.get(i);
				if(count>0){
					outstream.print(comma+i);
					comma=",";
				}
			}
		}

		if(counts_underscore!=null){
			outstream.println();
			outstream.println("Length_underscore counts:");

			for(int i=0; i<counts_underscore.length(); i++){
				int count=counts_underscore.get(i);
				if(count>0){outstream.println(i+"\t"+count);}
			}
		}

		if(counts_underscore2!=null){
			outstream.println();
			outstream.println("Length_underscore2 counts:");

			for(int i=0; i<counts_underscore2.length(); i++){
				int count=counts_underscore2.get(i);
				if(count>0){outstream.println(i+"\t"+count);}
			}
		}
		outstream.println();
		Shared.printMemory();
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
		
		LOADED=true;
	}
	
	/** Spawn process threads */
	private void spawnThreads(){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=Shared.threads();
		
		//Fill a list with ProcessThreads
		ArrayList<HashThread> alht=new ArrayList<HashThread>(ffin.length);
		for(int i=0; i<ffin.length; i++){
			alht.add(new HashThread(ffin[i]));
		}
		
		//Start the threads
		for(HashThread pt : alht){
			pt.start();
		}
		
		//Wait for completion of all threads
		boolean success=true;
		for(HashThread pt : alht){
			
			//Wait until this thread has terminated
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					//Attempt a join operation
					pt.join();
				} catch (InterruptedException e) {
					//Potentially handle this, if it is expected to occur
					e.printStackTrace();
				}
			}
			
			linesProcessed+=pt.linesProcessedT;
			linesValid+=pt.linesValidT;
			bytesProcessed+=pt.bytesProcessedT;
			success&=pt.success;
		}
		
		//Track whether any threads failed
		if(!success){errorState=true;}
		
		//Do anything necessary after processing
		System.gc();
	}
	
	/*--------------------------------------------------------------*/
	
	public static int get(String accession){
		if(accession==null){return -1;}
//		if(STRIP_UNDERSCORE){
//			accession=accession.replaceAll("[_-]", "");
//		}
		
		int dot=accession.indexOf('.');
		int len=(dot<0 ? accession.length() : dot);
		if(len>12 || !USE_TABLES){
			if(dot>-1){accession=accession.substring(0,  dot);}
			if(accession.length()<1){return -1;}
			int way=accession.charAt(0);
			Integer value=maps[way].get(accession);
			return value==null ? -1 : value.intValue();
		}else{
			long number=hash(accession);

			int value=tables.getCount(number);
			return value<1 ? -1 : value;
		}
	}
	
	private static long hash(String accession){
		long number=0;
		for(int i=0, max=accession.length(); i<max; i++){
			long c=accession.charAt(i);
			if(c=='.'){break;}
			if(c>='0' || c<='9'){c=c-'0';}
			else if(c=='_' || c=='-'){c=10;}//Collision, but should be OK
			else if(c>='A' && c<='Z'){c=c-'A'+11;}
			else{
				assert(false) : accession;
			}
			number=(number*37)+c;
		}
		return number;
	}
	
	public static int parseLineToTaxid(final byte[] line, final byte delimiter){
		int a=0, b=0;
		
		final int ncbi;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		assert(b>a) : "Missing field 0: "+new String(line);
		b++;
		a=b;
		
		while(b<line.length && line[b]!=delimiter){b++;}
//		assert(b>a) : "Missing field 1: "+new String(line);
		assert(b>=a) : "Missing field 1: "+new String(line)+"\n"+a+", "+b;
		//accession2=new String(line, a, b-a);
		b++;
		a=b;
		
		while(b<line.length && line[b]!=delimiter){b++;}
		assert(b>a) : "Missing field 2: "+new String(line);
		ncbi=Tools.parseInt(line, a, b);
		b++;
		a=b;
		
		return ncbi;
	}
	
	/*--------------------------------------------------------------*/
	
	public class HashThread extends Thread {
		
		@SuppressWarnings("unchecked")
		public HashThread(FileFormat ff_){
//			if(USE_MAPS){
				mapsT=new HashMap[128];
				for(int i=0; i<mapsT.length; i++){
					mapsT[i]=new HashMap<String, Integer>();
				}
//			}
			if(USE_TABLES){
				table=new HashBuffer(tables.tables(), 1000, 31, true, true);
			}
			ff=ff_;
		}
		
		@Override
		public void run(){
			
			ByteFile bf=ByteFile.makeByteFile(ff);
			
			byte[] line=bf.nextLine();
			while(line!=null && Tools.startsWith(line, "accession")){line=bf.nextLine();}
			
			while(line!=null){
				if(line.length>0){
					linesProcessedT++;
					bytesProcessedT+=line.length;
					
					final boolean valid=!Tools.startsWith(line, "accession");
					assert(valid);
					
					if(valid){
						boolean b=parseLine(line, (byte)'\t');
						if(b){linesValidT++;}
					}
				}
				line=bf.nextLine();
			}
			
			boolean closedError=bf.close();
			
//			if(USE_MAPS){
				for(int i=0; i<mapsT.length; i++){
					if(mapsT[i].size()>0){
						synchronized(maps[i]){
							maps[i].putAll(mapsT[i]);
						}
					}
					mapsT[i]=null;
				}
//			}
			if(USE_TABLES){
				long temp=table.flush();
			}
			
			success=!closedError;
		}
		
//		public boolean parseLineNumeric(final byte[] line, final byte delimiter){
//			int a=0, b=0;
//			
//			long accession=0;
//			final int ncbi, gi;
//			
//			while(b<line.length && line[b]!=delimiter){b++;}
//			assert(b>a) : "Missing field 0: "+new String(line);
//			for(int i=a; i<b; i++){
//				long c=line[i];
//				if(c=='.'){break;}
//				if(c<='9'){c=c-'0';}
//				else{c=c-'A'+10;}
//				accession=(accession*36)+c;
//			}
//			b++;
//			a=b;
//			
//			while(b<line.length && line[b]!=delimiter){b++;}
//			assert(b>a) : "Missing field 1: "+new String(line);
//			//accession2=new String(line, a, b-a);
//			b++;
//			a=b;
//			
//			while(b<line.length && line[b]!=delimiter){b++;}
//			assert(b>a) : "Missing field 2: "+new String(line);
//			ncbi=Tools.parseInt(line, a, b);
//			b++;
//			a=b;
//			
////			while(b<line.length && line[b]!=delimiter){b++;}
////			assert(b>a) : "Missing field 3: "+new String(line);
//////			gi=Tools.parseInt(line, a, b);
////			b++;
////			a=b;
//			
//			if(ncbi<1){return false;}
//			
//			if(tree!=null){
//				if(ncbi>=tree.nodes.length){return false;}
//				TaxNode tn=tree.getNode(ncbi);
//				if(tn==null || tn.level==TaxTree.NO_RANK || tn.level==TaxTree.LIFE || tn.level==TaxTree.DOMAIN){return false;}
//				if(tn.pid>=tree.nodes.length){return false;}
//				tn=tree.getNode(tn.pid);
//				if(tn==null || tn.level==TaxTree.NO_RANK || tn.level==TaxTree.LIFE){return false;}
//			}
//			assert(accession>=0) : new String(line);
//			table.set(accession, ncbi);
//			return true;
//		}
		
		public boolean parseLine(final byte[] line, final byte delimiter){
			int a=0, b=0;
			
			String accession;
			final int ncbi, gi;
			
			while(b<line.length && line[b]!=delimiter){b++;}
			assert(b>a) : "Missing field 0: "+new String(line);
			accession=new String(line, a, b-a);
			final int dot=accession.indexOf('.');
			if(dot>=0){//Should never happen
//				System.err.println(accession);
//				assert(dot==accession.length()-2) : accession;
				accession=accession.substring(0, dot);
			}
//			if(STRIP_UNDERSCORE){
//				accession=accession.replaceAll("[_-]", "");
//			}
			if(lengthCounts!=null){lengthCounts.incrementAndGet(b-a);}
			if(symbolCounts!=null){
				for(int i=a; i<b; i++){symbolCounts.incrementAndGet(line[i]);}
			}
			if(counts_underscore!=null && accession.contains("_")){counts_underscore.incrementAndGet(b-a);}
			if(counts_underscore2!=null && accession.indexOf('_')==2){counts_underscore2.incrementAndGet(b-a);}
			b++;
			a=b;
			
			while(b<line.length && line[b]!=delimiter){b++;}
//			assert(b>a) : "Missing field 1: "+new String(line);
			assert(b>=a) : "Missing field 1: "+new String(line)+"\n"+a+", "+b;
			//accession2=new String(line, a, b-a);
			b++;
			a=b;
			
			while(b<line.length && line[b]!=delimiter){b++;}
			assert(b>a) : "Missing field 2: "+new String(line);
			ncbi=Tools.parseInt(line, a, b);
			b++;
			a=b;
			
//			while(b<line.length && line[b]!=delimiter){b++;}
//			assert(b>a) : "Missing field 3: "+new String(line);
////			gi=Tools.parseInt(line, a, b);
//			b++;
//			a=b;
			
			if(ncbi<1){return false;}
			
			if(tree!=null){
				if(ncbi>=tree.nodes.length){return false;}
				TaxNode tn=tree.getNode(ncbi);
				if(tn==null || tn.levelExtended==TaxTree.NO_RANK_E || tn.levelExtended==TaxTree.LIFE_E || tn.levelExtended==TaxTree.DOMAIN_E){return false;}
				if(tn.pid>=tree.nodes.length){return false;}
				tn=tree.getNode(tn.pid);
				if(tn==null || tn.levelExtended==TaxTree.NO_RANK_E || tn.levelExtended==TaxTree.LIFE_E){return false;}
			}
			
			if(accession.length()<13 && USE_TABLES){
				long number=hash(accession);
				assert(number>=0) : new String(line);
				table.set(number, ncbi);
				return true;
			}
			
			int way=accession.charAt(0);
			mapsT[way].put(accession, ncbi);
//			Integer old=mapsT[way].put(accession, ncbi);
//			assert(old==null || old==ncbi) : "'"+accession+"': "+old+" -> "+ncbi;
//			System.err.println("'"+accession+"': "+old+" -> "+ncbi);
//			assert(dot==-1) : "'"+accession+"': "+old+" -> "+ncbi;
			return true;
		}
		
		private long linesProcessedT=0;
		private long linesValidT=0;
		private long bytesProcessedT=0;
		
		final FileFormat ff;
		HashMap<String, Integer>[] mapsT;
		HashBuffer table;
		boolean success=false;
		
	}
	
	/*--------------------------------------------------------------*/
	
	private void printOptions(){assert(false) : "Please read the associated shell script for usage information.";}
	
	
	/*--------------------------------------------------------------*/
	
	private ArrayList<String> in=new ArrayList<String>();
	private String out=null;
	
	/*--------------------------------------------------------------*/
	
	private long linesProcessed=0;
	private long linesValid=0;
	private long bytesProcessed=0;
	
	private AtomicIntegerArray lengthCounts=null;//new AtomicIntegerArray(20);
	private AtomicIntegerArray symbolCounts=null;//new AtomicIntegerArray(255);
	private AtomicIntegerArray counts_underscore=null;//new AtomicIntegerArray(20);
	private AtomicIntegerArray counts_underscore2=null;//new AtomicIntegerArray(20);
	
	/*--------------------------------------------------------------*/

	private final FileFormat ffin[];
	private final FileFormat ffout;
	
	
	/*--------------------------------------------------------------*/
	
	public static boolean LOADED=false;
	private static HashMap<String, Integer>[] maps=null;
	private static KmerTableSet tables;
	public static TaxTree tree=null;
//	public static final boolean USE_MAPS=true;
	public static final boolean USE_TABLES=true;
//	public static boolean STRIP_UNDERSCORE=false;
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=false;
	
}
