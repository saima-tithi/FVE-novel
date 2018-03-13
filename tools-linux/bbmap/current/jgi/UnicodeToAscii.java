package jgi;

import java.io.File;
import java.io.UnsupportedEncodingException;
import java.util.Arrays;

import fileIO.ReadWrite;
import fileIO.TextFile;
import fileIO.TextStreamWriter;
import shared.Parser;
import shared.ReadStats;
import shared.Shared;
import shared.Timer;
import shared.Tools;

/**
 * @author Brian Bushnell
 * @date Apr 21, 2015
 *
 */
public class UnicodeToAscii {
	
	public static void main(String[] args){
//		try {
//			System.err.println("AêñüC");
//			System.err.println(new String("AêñüC".getBytes("UTF8"), "UTF8"));
//		} catch (Exception e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
//		assert(false);
		
		Timer t=new Timer();
		t.start();
		UnicodeToAscii rr=new UnicodeToAscii(args);
		rr.process(t);
	}
	
	public UnicodeToAscii(String[] args){
		
		args=Parser.parseConfig(args);
		if(Parser.parseHelp(args, true)){
			assert(false) : "No help available.";
			System.exit(0);
		}
		
		System.err.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		
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
			
			if(parser.parse(arg, a, b)){
				//do nothing
			}else if(a.equals("null")){
				// do nothing
			}else if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
			}else if(parser.in1==null && i==0 && !arg.contains("=") && (arg.toLowerCase().startsWith("stdin") || new File(arg).exists())){
				parser.in1=arg;
			}else if(parser.out1==null && i==1 && !arg.contains("=")){
				parser.out1=arg;
			}else{
				System.err.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			in1=parser.in1;
			in2=parser.in2;

			out1=parser.out1;
			out2=parser.out2;
			
			overwrite=parser.overwrite;
			append=parser.append;
		}
		
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}
		if(out1!=null && out2==null && out1.indexOf('#')>-1){
			out2=out1.replace("#", "2");
			out1=out1.replace("#", "1");
		}
		
		if(in1==null){
			throw new RuntimeException("Error - at least one input file is required.");
		}
		
		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		if(out2!=null && out2.equalsIgnoreCase("null")){out2=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2)){
			System.err.println((out1==null)+", "+(out2==null)+", "+out1+", "+out2);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+"\n");
		}
		if(!Tools.testForDuplicateFiles(true, in1, in2, out1, out2) || !ReadStats.testFiles(false)){
			throw new RuntimeException("Duplicate filenames are not allowed.");
		}
	}
	
	private void process(Timer t){

		if(in1!=null && out1!=null){process(in1, out1);}
		if(in2!=null && out2!=null){process(in2, out2);}
		
	}
		
	private void process(String infile, String outfile){
		TextFile tf=new TextFile(infile, true);
		TextStreamWriter tsw=new TextStreamWriter(outfile, overwrite, append, true);
		tsw.start();
		for(String line=tf.readLine(false); line!=null; line=tf.readLine(false)){
			String line2=line;
			try {
				line2=new String(line.getBytes(), "UTF-8");
			} catch (UnsupportedEncodingException e) {
				try {
					line2=new String(line.getBytes(), "UTF-16");
				} catch (UnsupportedEncodingException e1) {}
			}
			tsw.println(Tools.fixHeader(line2));
//			tsw.println(Normalizer.normalize(line, Normalizer.Form.NFD));
		}
		tf.close();
		tsw.poisonAndWait();
	}

	private String in1, in2;
	private String out1, out2;
	@SuppressWarnings("unused")
	private boolean verbose=false;
	private boolean overwrite=true;
	private boolean append=false;
	
}
