package var2;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import shared.Parser;
import shared.Shared;
import shared.Tools;
import shared.Timer;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import stream.ConcurrentGenericReadInputStream;
import stream.FastaReadInputStream;

/**
 * @author Brian Bushnell
 * @date January 14, 2017
 *
 */
public class FilterVCF {
	
	public static void main(String[] args){
		Timer t=new Timer();
		FilterVCF sample=new FilterVCF(args);
		sample.process(t);
	}
	
	public FilterVCF(String[] args){
		
		args=Parser.parseConfig(args);
		if(Parser.parseHelp(args, true)){
			printOptions();
			System.exit(0);
		}
		
		outstream.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=Shared.threads();
		
		varFilter.clear();

		boolean setSamFilter=false;
		boolean setVarFilter=false;
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b==null || b.equalsIgnoreCase("null")){b=null;}
			while(a.startsWith("-")){a=a.substring(1);} //In case people use hyphens

			if(a.equals("lines")){
				maxLines=Long.parseLong(b);
				if(maxLines<0){maxLines=Long.MAX_VALUE;}
			}else if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("ref")){
				ref=b;
			}else if(a.equals("ploidy")){
				ploidy=Integer.parseInt(b);
			}else if(a.equals("sub")){
				Var.CALL_SUB=Tools.parseBoolean(b);
			}else if(a.equals("del")){
				Var.CALL_DEL=Tools.parseBoolean(b);
			}else if(a.equals("ins")){
				Var.CALL_INS=Tools.parseBoolean(b);
			}else if(a.equals("indel")){
				Var.CALL_INS=Var.CALL_DEL=Tools.parseBoolean(b);
			}
			
			else if(samFilter.parse(arg, a, b)){
				setSamFilter=true;
			}else if(varFilter.parse(a, b, arg)){
				setVarFilter=true;
			}
			
			else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			in1=parser.in1;
			out1=parser.out1;
			overwrite=parser.overwrite;
			append=parser.append;
		}

		if(!setSamFilter){samFilter=null;}
		if(!setVarFilter){varFilter=null;}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){
			printOptions();
			throw new RuntimeException("Error - at least two input files are required.");
		}
		
		if(!ByteFile.FORCE_MODE_BF2){
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF1=true;
		}
		

		
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out1+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1, ref)){
			throw new RuntimeException("\nCan't read to some input files.\n");
		}
		
		ffin1=FileFormat.testInput(in1, FileFormat.TXT, null, true, true);
		
		ffout1=FileFormat.testOutput(out1, FileFormat.TXT, null, true, overwrite, append, false);
		
		if(ref!=null){ScafMap.loadReference(ref, scafMap, samFilter);}
	}
	
	public void loadHeaderInScafMap(){
		for(byte[] line : header){
			if(Tools.startsWith(line, "##contig=<ID=")){
				scafMap.addFromVcf(line);
			}
		}
	}
	
	public String headerToString(){
		StringBuilder sb=new StringBuilder();
		for(byte[] line : header){
			for(byte b : line){
				sb.append((char)b);
			}
			sb.append('\n');
		}
		return sb.toString();
	}
	
	public void filter(FileFormat ff, ByteStreamWriter bsw){
		
		ByteFile bf=ByteFile.makeByteFile(ffin1);
		
		byte[] line=bf.nextLine();
		
		while(line!=null && (line.length==0 || line[0]=='#')){
			if(line.length>0){
				if(maxLines>0 && linesProcessed>=maxLines){break;}
				linesProcessed++;
				headerLinesProcessed++;
				bytesProcessed+=line.length;
				if(bsw!=null){bsw.println(line);}
				headerLinesOut++;
				header.add(line);
				if(Tools.startsWith(line, VCFFile.CHROM_POS)){
					String[] split=new String(line).split("\t");
					for(int i=9; i<split.length; i++){
						samples.add(split[i]);
					}
				}else{
					String[] split=new String(line).split("=");
					if(split.length==2){
						String a=split[0], b=split[1];
						if(a.equalsIgnoreCase("##ploidy")){
							ploidy=Integer.parseInt(b);
						}else if(a.equalsIgnoreCase("##pairingRate")){
							properPairRate=(float) Double.parseDouble(b);
						}else if(a.equalsIgnoreCase("##totalQualityAvg")){
							totalQualityAvg=(float) Double.parseDouble(b);
						}else if(a.equalsIgnoreCase("##mapqAvg")){
							totalMapqAvg=(float) Double.parseDouble(b);
						}else if(a.equalsIgnoreCase("#readLengthAvg")){
							readLengthAvg=(float) Double.parseDouble(b);
						}
					}
				}
			}
			line=bf.nextLine();
		}
		
		loadHeaderInScafMap();
		if(ScafMap.defaultScafMap==null){ScafMap.defaultScafMap=scafMap;}
		assert(ScafMap.defaultScafMap.size()>0) : ScafMap.defaultScafMap+"\n"+headerToString();
		boolean varFormatOK=true;
		while(line!=null){
			if(line.length>0){
				if(maxLines>0 && linesProcessed>=maxLines){break;}
				linesProcessed++;
				bytesProcessed+=line.length;
				
				final boolean isHeader=(line[0]=='#');
				
				if(isHeader){
					assert(false) : "Encountered intermediate header.";
					headerLinesProcessed++;
					if(bsw!=null){bsw.println(line);}
					header.add(line);
					if(Tools.startsWith(line, VCFFile.CHROM_POS)){
						String[] split=new String(line).split("\t");
						for(int i=9; i<split.length; i++){
							samples.add(split[i]);
						}
					}
				}else{
					variantLinesProcessed++;
					VCFLine vline=new VCFLine(line);
					boolean pass=true;
					if(pass && samFilter!=null){pass&=samFilter.passesFilter(vline);}
					if(pass && varFilter!=null){
						Var v=null;
						
						if(varFormatOK){
							try {
								v=vline.toVar();
							} catch (Throwable e) {
								System.err.println("WARNING: This VCF file does not support Var format.\n"
										+ "Filtering can only be done on location and quality score.\n");
								varFormatOK=false;
							}
						}

						if(!Var.CALL_DEL && vline.type()==Var.DEL){pass=false;}
						else if(!Var.CALL_INS && vline.type()==Var.INS){pass=false;}
						else if(!Var.CALL_SUB && vline.type()==Var.SUB){pass=false;}
						
						if(v!=null){
							pass&=varFilter.passesFilter(v, properPairRate, 
									totalQualityAvg, totalMapqAvg, readLengthAvg, ploidy, scafMap);
						}else{
							pass&=vline.qual>=varFilter.minScore;
						}
					}
					if(pass){
						if(bsw!=null){bsw.println(line);}
						variantLinesOut++;
					}
				}
			}
			line=bf.nextLine();
		}

		errorState|=bf.close();
		if(bsw!=null){errorState|=bsw.poisonAndWait();}
	}
	
	void process(Timer t){
		
		ByteStreamWriter bsw;
		if(ffout1!=null){
			bsw=new ByteStreamWriter(ffout1);
			bsw.start();
		}else{bsw=null;}
		
		filter(ffin1, bsw);
		
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
		outstream.println("Header Lines In:   \t"+headerLinesProcessed);
		outstream.println("Variant Lines In:  \t"+variantLinesProcessed);
		outstream.println("Header Lines Out:  \t"+headerLinesOut);
		outstream.println("Variant Lines Out: \t"+variantLinesOut);
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	
	private void printOptions(){assert(false) : "Please read the associated shell script for usage information.";}
	
	/*--------------------------------------------------------------*/

	private long linesProcessed=0;
	private long headerLinesProcessed=0;
	private long variantLinesProcessed=0;
	private long headerLinesOut=0;
	private long variantLinesOut=0;
	private long bytesProcessed=0;
	
	private long maxLines=Long.MAX_VALUE;

	public ArrayList<byte[]> header=new ArrayList<byte[]>();
	public ArrayList<String> samples=new ArrayList<String>();
	
	SamFilter samFilter=new SamFilter();
	VarFilter varFilter=new VarFilter();
	
	/*--------------------------------------------------------------*/
	
	public int ploidy=1;
	public float properPairRate=0;
	public float totalQualityAvg=30;
	public float totalMapqAvg=30;
	public float readLengthAvg=150;
	
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	private String out1=null;
	private String ref=null;

	private final FileFormat ffin1;
	private final FileFormat ffout1;
	
	public final ScafMap scafMap=new ScafMap();
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=false;
	private boolean append=false;
	
}
