package tax;

import java.io.File;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.Locale;

import shared.Parser;
import shared.ReadStats;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ByteBuilder;
import stream.FastaReadInputStream;
import dna.Data;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;

/**
 * @author Brian Bushnell
 * @date April 4, 2017
 *
 */
public class ShrinkAccession {
	
	public static void main(String[] args){
		Timer t=new Timer();
		ShrinkAccession mb=new ShrinkAccession(args);
		mb.process(t);
	}
	
	public ShrinkAccession(String[] args){
		
		args=Parser.parseConfig(args);
		if(Parser.parseHelp(args, true)){
			printOptions();
			System.exit(0);
		}
		
		outstream.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=Shared.threads();
		if(Data.PIGZ()){
			ReadWrite.ZIPLEVEL=Tools.max(ReadWrite.ZIPLEVEL, 6);
		}
		
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
			}else if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				ReadWrite.verbose=verbose;
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
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			
			in1=parser.in1;

			out1=parser.out1;
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

		ffout1=FileFormat.testOutput(out1, FileFormat.TXT, null, true, overwrite, append, false);
		ffin1=FileFormat.testInput(in1, FileFormat.TXT, null, true, true);
		
	}
	
	void process(Timer t){
		
		ByteFile bf=ByteFile.makeByteFile(ffin1);
		ByteStreamWriter bsw=new ByteStreamWriter(ffout1);
		bsw.start();

		long linesProcessed=0;
		long charsProcessed=0;
		long badLines=0;
		
		byte[] line=bf.nextLine();
		boolean valid=false;
		ByteBuilder bb=new ByteBuilder(10000);
		while(line!=null){
			if(Tools.startsWith(line, "accession\t")){
				bb.append(line);
				bb.append('\n');
			}else{
				charsProcessed+=line.length+1;
				linesProcessed++;
				
				if(AccessionToTaxid.parseLineToTaxid(line, (byte)'\t')<1){
					badLines++;
				}else{
					int i=0;
					while(i<line.length){
						byte b=line[i];
						bb.append(b);
						i++;
						if(b=='\t'){break;}
					}
					while(i<line.length){
						byte b=line[i];
//						bb.append(b);
						i++;
						if(b=='\t'){
							bb.append(b);
							break;
						}
					}
					while(i<line.length){
						byte b=line[i];
						bb.append(b);
						i++;
						if(b=='\t'){break;}
					}
					bb.append('\n');
				}
				
//				String[] split=new String(line).split("\t");
//				bb.append(split[0]);
//				bb.append('\t');
//				bb.append('\t');
//				bb.append(split[2]);
//				bb.append('\t');
//				bb.append('\n');
			}
			if(bb.length()>8000){
				bsw.print(bb);
				bb.clear();
			}
			line=bf.nextLine();
		}
		if(bb.length()>0){
			bsw.print(bb);
			bb.clear();
		}
		
		errorState|=bf.close();
		if(bsw!=null){errorState|=bsw.poisonAndWait();}
		
		t.stop();
		
		double rpnano=linesProcessed/(double)(t.elapsed);
		double bpnano=charsProcessed/(double)(t.elapsed);

		String rpstring=(linesProcessed<100000 ? ""+linesProcessed : linesProcessed<100000000 ? (linesProcessed/1000)+"k" : (linesProcessed/1000000)+"m");
		String bpstring=(charsProcessed<100000 ? ""+charsProcessed : charsProcessed<100000000 ? (charsProcessed/1000)+"k" : (charsProcessed/1000000)+"m");

		while(rpstring.length()<8){rpstring=" "+rpstring;}
		while(bpstring.length()<8){bpstring=" "+bpstring;}
		
		outstream.println("Discarded "+badLines+" lines.\n");
		
		outstream.println("Time:                         \t"+t);
		outstream.println("Reads Processed:    "+rpstring+" \t"+String.format(Locale.ROOT, "%.2fk lines/sec", rpnano*1000000));
		outstream.println("Bases Processed:    "+bpstring+" \t"+String.format(Locale.ROOT, "%.2fm chars/sec", bpnano*1000));
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	
	private void printOptions(){assert(false) : "Please read the associated shell script for usage information.";}
	
	
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	private String out1=null;
	
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
