package sketch;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import dna.AminoAcid;
import stream.ByteBuilder;
import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import structures.ListNum;
import structures.LongHashSet;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parser;
import shared.ReadStats;
import shared.Shared;
import shared.Timer;
import shared.Tools;

/**
 * @author Brian Bushnell
 * @date Oct 17, 2014
 *
 */
public class InvertKey {
	
	public static void main(String[] args){
		Timer t=new Timer();
		InvertKey mb=new InvertKey(args);
		mb.process(t);
	}
	
	public InvertKey(String[] args){
		
		args=Parser.parseConfig(args);
		if(Parser.parseHelp(args, true)){
			printOptions();
			System.exit(0);
		}
		
		for(String s : args){if(s.startsWith("out=standardout") || s.startsWith("out=stdout")){outstream=System.err;}}
		outstream.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");

		
		
		Shared.capBufferLen(200);
		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=Shared.threads();
		
		FASTQ.TEST_INTERLEAVED=FASTQ.FORCE_INTERLEAVED=false;
		int k_=31, k2_=0;
		
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
//				align2.FastaReadInputStream2.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("key")){
				keyString=b;
			}else if(a.equals("out")){
				out1=b;
			}else if(a.equalsIgnoreCase("k")){
				if(b.indexOf(',')>=0){
					String[] bsplit=b.split(",");
					assert(bsplit.length==2) : "Bad argument "+arg;
					int x=Integer.parseInt(bsplit[0]);
					int y=Integer.parseInt(bsplit[1]);
					k_=Tools.max(x, y);
					k2_=Tools.min(x, y);
				}else{
					k_=Integer.parseInt(b);
					k2_=0;
				}
			}else if(a.equalsIgnoreCase("k2")){
				k2_=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("printonce")){
				printOnce=Tools.parseBoolean(b);
			}else if(parser.in1==null && i==0 && !arg.contains("=") && (arg.toLowerCase().startsWith("stdin") || new File(arg).exists())){
				parser.in1=arg;
			}else if(parser.out1==null && i==1 && !arg.contains("=")){
				out1=arg;
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		SketchObject.k=k=k_;
		SketchObject.k2=k2=k2_;
		shift=2*k;
		shift2=shift-2;
		mask=(shift>63 ? -1L : ~((-1L)<<shift)); //Conditional allows K=32
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			
			in1=parser.in1;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){
			printOptions();
			throw new RuntimeException("Error - at least one input file is required.");
		}
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2){
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF1=true;
		}

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+"\n");
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTA, null, true, overwrite, append, false);

		ffin1=FileFormat.testInput(in1, FileFormat.FASTA, null, true, true);
		
		SketchObject.postParse();
		
		if(keyString.indexOf(',')>0){
			String[] split=keyString.split(",");
			set=new LongHashSet(split.length*2);
			for(String s : split){
				long x=Long.MAX_VALUE-Sketch.parseA48(s);
				set.add(x);
//				assert(set.contains(x)) : x+", "+set.size()+", "+set.toStringListView();
			}
			key0=-1;
//			System.err.println(set.toStringListView()+", "+set.size());
			assert(!set.isEmpty());
		}else if(keyString.endsWith(".sketch")){
			SketchTool tool=new SketchTool(10000, 0, false);
			Sketch sk=tool.loadSketches(keyString, null, SketchObject.ONE_SKETCH, 1, 1000000, 0).get(0);
			set=new LongHashSet(sk.length()*2);
			for(long x : sk.array){set.add(Long.MAX_VALUE-x);}
			key0=-1;
//			System.err.println(set.toStringListView()+", "+set.size());
			assert(!set.isEmpty());
		}else{
			key0=Long.MAX_VALUE-Sketch.parseA48(keyString);
			set=null;
//			System.err.println(key0);
		}
	}
	
	void process(Timer t){
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, null, null, null);
			cris.start();
			if(verbose){outstream.println("Started cris");}
		}
		boolean paired=cris.paired();
//		if(verbose){
			if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
//		}

		final ByteStreamWriter bsw;
		if(out1!=null){
			fasta=ffout1.fasta() && !out1.endsWith(".txt");
			bsw=new ByteStreamWriter(ffout1);
			bsw.start();
		}else{bsw=null;}
		
		long readsProcessed=0;
		long basesProcessed=0;
		boolean finished=false;
		
		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
//			outstream.println("Fetched "+reads);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}
			
			while(reads!=null && reads.size()>0 && !finished){
				
				for(int idx=0; idx<reads.size() && !finished; idx++){
					final Read r1=reads.get(idx);

					finished=invert(key0, r1, bsw);
					
					final int initialLength1=r1.length();
					
					readsProcessed++;
					basesProcessed+=initialLength1;
				}

				cris.returnList(ln.id, ln.list.isEmpty());
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		errorState|=(ReadWrite.closeStream(cris));
		if(bsw!=null){errorState|=bsw.poisonAndWait();}
		
		t.stop();
		
		double rpnano=readsProcessed/(double)(t.elapsed);
		double bpnano=basesProcessed/(double)(t.elapsed);

		String rpstring=(readsProcessed<100000 ? ""+readsProcessed : readsProcessed<100000000 ? (readsProcessed/1000)+"k" : (readsProcessed/1000000)+"m");
		String bpstring=(basesProcessed<100000 ? ""+basesProcessed : basesProcessed<100000000 ? (basesProcessed/1000)+"k" : (basesProcessed/1000000)+"m");

		while(rpstring.length()<8){rpstring=" "+rpstring;}
		while(bpstring.length()<8){bpstring=" "+bpstring;}
		
		outstream.println("Time:                         \t"+t);
		outstream.println("Reads Processed:    "+rpstring+" \t"+String.format(Locale.ROOT, "%.2fk reads/sec", rpnano*1000000));
		outstream.println("Bases Processed:    "+bpstring+" \t"+String.format(Locale.ROOT, "%.2fm bases/sec", bpnano*1000));
		
		if(errorState && !finished && maxReads<1){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	private boolean invert(long key2, Read r, ByteStreamWriter bsw) {
		final byte[] bases=r.bases;
		
		long kmer=0;
		long rkmer=0;
		int len=0;
		

//		System.err.println("Looking for "+key+"\t"+Sketch.toA48(key)+"\t"+Sketch.toA48(Long.MAX_VALUE-key));
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			long x=AminoAcid.baseToNumber[b];
			long x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=(rkmer>>>2)|(x2<<shift2);
			if(x<0){len=0;}else{len++;}
			if(len>=k){
				kmersProcessed++;
				long z=Tools.max(kmer, rkmer);
				long key=SketchTool.hash(z);
				boolean found=(key0>=0 ? key==key0 : set.contains(key));
				if(found){
					if(fasta){bsw.println(">"+Sketch.toA48(Long.MAX_VALUE-key)+" "+(i-k+1)+" "+r.id);}
					bsw.println(AminoAcid.kmerToString(Tools.min(kmer, rkmer), k));
					if(printOnce){
						if(key0>=0){return true;}
						else{
							set.remove(key);
							return set.isEmpty();
						}
					}
				}
			}
		}
		return false;
	}

	/*--------------------------------------------------------------*/
	
	private void printOptions(){
		throw new RuntimeException("Please run the shellscript with no arguments for usage information."); //TODO
	}
	
	
	/*--------------------------------------------------------------*/
	
	final long key0;
	final LongHashSet set;
	
	final int shift;
	final int shift2;
	final int k;
	final int k2;
	final long mask;
	
	boolean printOnce=true;
	long kmersProcessed=0;
	
	private String in1=null;
	boolean fasta;
	boolean sketch;
	private String keyString=null;

	private String out1="stdout.fa";
	
	/*--------------------------------------------------------------*/

	private long maxReads=-1;
	
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffin1;

	private final FileFormat ffout1;
	
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=false;
	private boolean append=false;
	
}
