package driver;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;
import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import shared.Parser;
import shared.ReadStats;
import shared.Shared;
import shared.Timer;
import shared.Tools;

/**
 * @author Brian Bushnell
 * @date Feb 21, 2017
 *
 */
public class PlotGC {
	
	public static void main(String[] args){
		Timer t=new Timer();
		PlotGC mb=new PlotGC(args);
		mb.process(t);
	}
	
	public PlotGC(String[] args){
		
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
				ReadWrite.verbose=verbose;
			}else if(a.equals("interval")){
				interval=Integer.parseInt(b);
			}else if(a.equals("offset")){
				offset=Integer.parseInt(b);
			}else if(a.equals("printshortbins") || a.equals("psb")){
				printShortBins=Tools.parseBoolean(b);
			}else if(a.equals("out")){
				out1=b;
			}else if(parser.parse(arg, a, b)){
				//do nothing 
			}else if(parser.in1==null && i==0 && !arg.contains("=") && (arg.toLowerCase().startsWith("stdin") || new File(arg).exists())){
				parser.in1=arg;
			}else if(parser.out1==null && i==1 && !arg.contains("=")){
				parser.out1=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			
			in1=parser.in1;
			
			extin=parser.extin;
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
		
		ffout1=FileFormat.testOutput(out1, FileFormat.TEXT, ".txt", true, overwrite, append, false);

		ffin1=FileFormat.testInput(in1, FileFormat.FASTA, extin, true, true);
	}
	
	void process(Timer t){
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, null, null, null);
			cris.start();
			if(verbose){outstream.println("Started cris");}
		}

		final TextStreamWriter tsw;
		if(out1!=null){
			tsw=new TextStreamWriter(ffout1);
			tsw.start();
			tsw.println("name\tinterval\tstart\tstop\trunningStart\trunningStop\tgc");
		}else{tsw=null;}
		
		long readsProcessed=0;
		long basesProcessed=0;
		
		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
//			outstream.println("Fetched "+reads);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}
			
			long rStart=0, rStop=0;
			
			int[] acgt=new int[4];
			
			while(reads!=null && reads.size()>0){
				
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					
					final int initialLength1=r1.length();
					
					readsProcessed++;
					basesProcessed+=initialLength1;
				}
				
				for(Read r : reads){
					Arrays.fill(acgt, 0);
					byte[] bases=r.bases;
					int next=interval-1;
					int start=0, i=0;
					for(; i<bases.length; i++){
						int num=AminoAcid.baseToNumber[bases[i]];
						if(num>=0){acgt[num]++;}
						if(i>=next){
							int len=(i-start);
							rStop=rStart+len;
							String s=toGC(r.id, start, i, rStart, rStop, acgt);
							start=i+1;
							rStart=rStop+1;
							next=i+interval;
							Arrays.fill(acgt, 0);
							if(tsw!=null && s!=null){
								tsw.print(s);
							}
						}
					}
					if(printShortBins && i>start){
						int len=(i-start)-1;
						rStop=rStart+len;
						String s=toGC(r.id, start, i-1, rStart, rStop, acgt);
						start=i+1;
						rStart=rStop+1;
						next=i+interval;
						Arrays.fill(acgt, 0);
						if(tsw!=null && s!=null){
							tsw.print(s);
						}
					}
				}
				
				cris.returnList(ln.id, ln.list.isEmpty());
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		errorState|=ReadWrite.closeStreams(cris);
		
		if(tsw!=null){
			errorState|=tsw.poisonAndWait();
		}
		
		t.stop();
		
		if(showSpeed){
			double rpnano=readsProcessed/(double)(t.elapsed);
			double bpnano=basesProcessed/(double)(t.elapsed);
			
			String rpstring=(readsProcessed<100000 ? ""+readsProcessed : readsProcessed<100000000 ? (readsProcessed/1000)+"k" : (readsProcessed/1000000)+"m");
			String bpstring=(basesProcessed<100000 ? ""+basesProcessed : basesProcessed<100000000 ? (basesProcessed/1000)+"k" : (basesProcessed/1000000)+"m");
			
			while(rpstring.length()<8){rpstring=" "+rpstring;}
			while(bpstring.length()<8){bpstring=" "+bpstring;}
			
			outstream.println("Time:                         \t"+t);
			outstream.println("Reads Processed:    "+rpstring+" \t"+String.format(Locale.ROOT, "%.2fk reads/sec", rpnano*1000000));
			outstream.println("Bases Processed:    "+bpstring+" \t"+String.format(Locale.ROOT, "%.2fm bases/sec", bpnano*1000));
		}
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	private String toGC(String name, int start, int stop, long rstart, long rstop, int[] acgt){
		int at=acgt[0]+acgt[3];
		int gc=acgt[1]+acgt[2];
		float sum=Tools.max(1, at+gc);
		float gcf=gc/sum;
		return String.format(Locale.ROOT, "%s\t%d\t%d\t%d\t%d\t%d\t%.3f\n", name, stop-start+1, start+offset, stop+offset, rstart+offset, rstop+offset, gcf);
	}
	
	/*--------------------------------------------------------------*/
	
	private void printOptions(){
		assert(false) : "printOptions: TODO";
	}
	
	
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	
	private String out1="stdout.txt";
	
	private String extin=null;
	
	/*--------------------------------------------------------------*/

	private long maxReads=-1;
	
	private int interval=1000;
	private int offset=0;

	private boolean showSpeed=false;
	private boolean printShortBins=true;
	
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
