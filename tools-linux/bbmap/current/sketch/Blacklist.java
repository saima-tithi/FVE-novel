package sketch;

import java.io.File;
import java.util.ArrayList;

import dna.Data;
import kmer.AbstractKmerTable;

public class Blacklist extends SketchObject {
	
	public static void addFiles(String fname){
		if(fname==null){return;}
		ArrayList<Sketch> sketches=new ArrayList<Sketch>();
		if(fname.indexOf(',')<0 || new File(fname).exists()){
			ArrayList<Sketch> temp=addFile(fname);
			sketches.addAll(temp);
		}else{
			String[] split=fname.split(",");
			for(String s : split){
				ArrayList<Sketch> temp=addFile(s);
				sketches.addAll(temp);
			}
		}
		addSketches(sketches);
	}
	
	private static ArrayList<Sketch> addFile(String fname){
		System.err.println("Adding "+fname+" to blacklist.");
		assert(!added.contains(fname));
		added.add(fname);
		SketchTool tool=new SketchTool(1000000, 1, false);
		ArrayList<Sketch> sketches=tool.loadSketches(fname, null, ONE_SKETCH, 1f, -1, defaultParams.minEntropy);
		return sketches;
	}
	
	private static void addSketches(ArrayList<Sketch> sketches){
		long size=0;
		for(Sketch sk : sketches){
			size+=sk.length();
		}
		long size2=(size*4)/3;
		assert(size2>0 && size2+1000<Integer.MAX_VALUE) : size2;
		if(keySets==null){
			keySets=AbstractKmerTable.preallocate(ways, AbstractKmerTable.ARRAY1D, (int)size2, true);
		}
		for(Sketch sk : sketches){
			for(long key : sk.array){
				increment(Long.MAX_VALUE-key);
			}
		}
	}
	
	public static int increment(long key){
		int way=0;//(int)(key%ways);
		return keySets[way].increment(key);
	}
	
	public static boolean contains(long key){
		if(keySets==null){return false;}
		int way=0;//(int)(key%ways);
		return keySets[way].getValue(key)>0;
	}
	
	public static boolean exists(){
		return keySets!=null;
	}

	static synchronized String ntBlacklist(){return ntBlacklist!=null ? ntBlacklist : (ntBlacklist=Data.findPath("?blacklist_nt_species_1000.sketch"));}
	static synchronized String silvaBlacklist(){return silvaBlacklist!=null ? silvaBlacklist : (silvaBlacklist=Data.findPath("?blacklist_silva_species_500.sketch"));}
	static synchronized String refseqBlacklist(){return refseqBlacklist!=null ? refseqBlacklist : (refseqBlacklist=Data.findPath("?blacklist_refseq_species_300.sketch"));}
	static synchronized String imgBlacklist(){return imgBlacklist!=null ? imgBlacklist : (imgBlacklist=Data.findPath("?blacklist_img_species_300.sketch"));}
	static synchronized String nrBlacklist(){return null;}//Data.findPath("?blacklist_nr_species_1000.sketch");

	private static String ntBlacklist;
	private static String silvaBlacklist;
	private static String refseqBlacklist;
	private static String imgBlacklist;
	private static String nrBlacklist;
	
	/** Hold kmers.  A kmer X such that X%WAYS=Y will be stored in keySets[Y] */
	private static AbstractKmerTable[] keySets;
	private static final int ways=1;
//	private static final int initialSize=16000;
	private static ArrayList<String> added=new ArrayList<String>();
	
}
