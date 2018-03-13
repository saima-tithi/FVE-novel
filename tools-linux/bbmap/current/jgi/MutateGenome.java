package jgi;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import java.util.Random;

import stream.ByteBuilder;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.Read;
import structures.ListNum;
import dna.AminoAcid;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parser;
import shared.Shared;
import shared.Timer;

/**
 * @author Brian Bushnell
 * @date June 1, 2016
 *
 */
public class MutateGenome {

	public static void main(String[] args){
		Timer t=new Timer();
		MutateGenome mg=new MutateGenome(args);
		mg.process(t);
	}
	
	public MutateGenome(String[] args){
		
		args=Parser.parseConfig(args);
		if(Parser.parseHelp(args, true)){
			printOptions();
			System.exit(0);
		}
		
		outstream.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		
		Shared.setBufferLen(1);
		FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b==null || b.equalsIgnoreCase("null")){b=null;}
			while(a.startsWith("-")){a=a.substring(1);} //In case people use hyphens

			if(parser.parse(arg, a, b)){
				//do nothing
			}else if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else if(a.equals("period")){
				period=Integer.parseInt(b);
			}else if(a.equals("subrate") || a.equals("snprate")){
				subRate=Float.parseFloat(b);
				if(subRate>1){subRate/=100;}
			}else if(a.equals("indelrate")){
				indelRate=Float.parseFloat(b);
				if(indelRate>1){indelRate/=100;}
			}else if(a.equals("seed")){
				seed=Long.parseLong(b);
			}else if(a.equals("prefix")){
				prefix=b;
			}else if(a.equals("id") || a.equals("identity")){
				float x=Float.parseFloat(b);
				if(x>1){x=x/100;}
				x=1-x;
				indelRate=x*0.01f;
				subRate=x*0.99f;
			}else if(a.equals("fraction") || a.equals("completeness")){
				float x=Float.parseFloat(b);
				if(x>1){x=x/100;}
				genomeFraction=x;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		errorRate=subRate+indelRate;

		assert(subRate>=0 && subRate<=1) : "Substitution rate must be between 0 and 1, inclusive.  Invalid value: "+subRate;
		assert(indelRate>=0 && indelRate<=1) : "Indel rate must be between 0 and 1, inclusive.  Invalid value: "+indelRate;
		assert(errorRate>=0 && errorRate<=1) : "Total error rate must be between 0 and 1, inclusive.  Invalid value: "+errorRate;
		
		System.err.println(String.format(Locale.ROOT, "Target Identity:   \t%.2f%%", (1-errorRate)*100));
		System.err.println(String.format(Locale.ROOT, "Substitution Rate: \t%.2f%%", subRate*100));
		System.err.println(String.format(Locale.ROOT, "Indel Rate:        \t%.2f%%", indelRate*100));
		
		if(seed<0){randy=new Random();}
		else{randy=new Random(seed);}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			in1=parser.in1;
			out1=parser.out1;
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTA, null, true, true, false, false);
		ffin1=FileFormat.testInput(in1, FileFormat.FASTA, null, true, true);
	}
	
	void process(Timer t){
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, null);
			cris.start();
		}

		final ConcurrentReadOutputStream ros;
		if(out1!=null){
			final int buff=4;
			
			if(cris.paired() && (in1==null || !in1.contains(".sam"))){
				outstream.println("Writing interleaved.");
			}			

			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1)) : "Input file and output file have same name.";
			
			ros=ConcurrentReadOutputStream.getStream(ffout1, null, buff, null, false);
			ros.start();
		}else{ros=null;}

		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}
			ByteBuilder bb=new ByteBuilder();

			while(reads!=null && reads.size()>0){
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}
				
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					
					readsProcessed++;
					basesProcessed+=r1.length();
					
					processRead(r1, bb);
				}
				
				if(ros!=null){ros.add(reads, ln.id);}

				cris.returnList(ln.id, ln.list.isEmpty());
				if(verbose){outstream.println("Returned a list.");}
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		ReadWrite.closeStreams(cris, ros);
		if(verbose){outstream.println("Finished.");}
		
		
		{
			t.stop();
			
			//Calculate units per nanosecond
			double rpnano=readsProcessed/(double)(t.elapsed);
			double bpnano=basesProcessed/(double)(t.elapsed);
			double brnano=basesRetained/(double)(t.elapsed);
			
			//Add "k" and "m" for large numbers
			String rpstring=(readsProcessed<100000 ? ""+readsProcessed : readsProcessed<100000000 ? (readsProcessed/1000)+"k" : (readsProcessed/1000000)+"m");
			String bpstring=(basesProcessed<100000 ? ""+basesProcessed : basesProcessed<100000000 ? (basesProcessed/1000)+"k" : (basesProcessed/1000000)+"m");
			String brstring=(basesRetained<100000 ? ""+basesRetained : basesRetained<100000000 ? (basesRetained/1000)+"k" : (basesRetained/1000000)+"m");
			String mastring=(mutationsAdded<100000 ? ""+mutationsAdded : mutationsAdded<100000000 ? (mutationsAdded/1000)+"k" : (mutationsAdded/1000000)+"m");
			
			//Format the strings so they have they are right-justified
			while(rpstring.length()<8){rpstring=" "+rpstring;}
			while(bpstring.length()<8){bpstring=" "+bpstring;}
			while(brstring.length()<8){brstring=" "+brstring;}
			while(mastring.length()<8){mastring=" "+mastring;}
			
			outstream.println("\nTime:                         \t"+t);
			outstream.println("Reads Processed:    "+rpstring+" \t"+String.format(Locale.ROOT, "%.2fk reads/sec", rpnano*1000000));
			outstream.println("Bases Processed:    "+bpstring+" \t"+String.format(Locale.ROOT, "%.2fm bases/sec", bpnano*1000));
			if(genomeFraction<1){outstream.println("Bases Retained:     "+brstring+" \t"+String.format(Locale.ROOT, "%.2fm bases/sec", brnano*1000));}
			outstream.println("Mutations Added:    "+mastring+" \t"+String.format(Locale.ROOT, "%.2f%% Identity", 100f-mutationsAdded*100f/basesProcessed));
		}
		
		t.stop();
	}
	
	public void processRead(Read r, ByteBuilder bb){
		
		//Setup
		bb.clear();
		r.quality=null;
		
		//Handle genomeFraction
		if(genomeFraction<1){
			final byte[] bases0=r.bases;
			int retain=(int)(bases0.length*(genomeFraction));
			if(retain<bases0.length){
				final int start=randy.nextInt(bases0.length);
				int i=0, j=start;
				for(; i<retain && j<bases0.length; i++, j++){
					bb.append(bases0[j]);
				}
				j=0;
				
				if(i<retain){mutationsAdded++;} //Chimeric junction
				
				for(; i<retain; i++, j++){
					bb.append(bases0[j]);
				}
				r.bases=bb.toBytes();
				bb.clear();
			}
		}
		
		//Handle mutations
		final byte[] bases=r.bases;
		
		if(period>-1){
			int basesSinceMutation=0;
			for(int i=0; i<bases.length; i++){
				byte b=bases[i];
				if(basesSinceMutation>=period && AminoAcid.isFullyDefined(b)){
					basesSinceMutation=0;
					mutationsAdded++;
					float x=randy.nextFloat()*errorRate;
					if(x<subRate){
						b=AminoAcid.numberToBase[((AminoAcid.baseToNumber[b]+randy.nextInt(3)+1)&3)];
						bb.append(b);
					}else if(randy.nextBoolean()){//del
						//do nothing
					}else{//ins
						i--;
						b=AminoAcid.numberToBase[randy.nextInt(4)];
						bb.append(b);
					}
				}else{
					basesSinceMutation++;
					bb.append(b);
				}
			}
		}else{
			for(int i=0; i<bases.length; i++){
				byte b=bases[i];
				float x=randy.nextFloat();
				if(x<errorRate && AminoAcid.isFullyDefined(b)){
					mutationsAdded++;
					if(x<subRate){
						b=AminoAcid.numberToBase[((AminoAcid.baseToNumber[b]+randy.nextInt(3)+1)&3)];
						bb.append(b);
					}else if(randy.nextBoolean()){//del
						//do nothing
					}else{//ins
						i--;
						b=AminoAcid.numberToBase[randy.nextInt(4)];
						bb.append(b);
					}
				}else{
					bb.append(b);
				}
			}
		}
		
		//Modify read
		r.bases=bb.toBytes();
		
		if(prefix!=null){
			r.id=prefix+r.numericID;
		}
		basesRetained+=r.bases.length;
	}
	
	/*--------------------------------------------------------------*/
	
	/** This is called if the program runs with no parameters */
	private void printOptions(){assert(false) : "Please read the associated shell script for usage information.";}
	
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	private String out1=null;

	private String prefix=null;
	
	private final FileFormat ffin1;
	private final FileFormat ffout1;
	
	/*--------------------------------------------------------------*/

	private long maxReads=-1;
	private long mutationsAdded=0;

	private int period=-1;
	
	private float genomeFraction=1;
	private long basesRetained;

	private long readsProcessed=0;
	private long basesProcessed=0;

	private float subRate=0;
	private float indelRate=0;
	private final float errorRate;
	
	private final Random randy;
	private long seed=-1;
	
	/*--------------------------------------------------------------*/
	
	private java.io.PrintStream outstream=System.err;
	public static boolean verbose=false;
	
}
