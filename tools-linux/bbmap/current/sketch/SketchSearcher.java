package sketch;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;

import shared.Parser;
import shared.Shared;
import shared.Tools;
import structures.AbstractBitSet;
import structures.Heap;
import tax.TaxNode;

public class SketchSearcher extends SketchObject {
	
	public SketchSearcher(){
		
	}

	public boolean parse(String arg, String a, String b, boolean addFileIfNotFound){
		if(Parser.isJavaFlag(arg)){
			//do nothing
		}else if(parseSketchFlags(arg, a, b)){
			//Do nothing
		}else if(defaultParams.parse(arg, a, b)){
			//Do nothing
		}else if(a.equals("verbose")){
			verbose=Tools.parseBoolean(b);
		}else if(a.equals("ref")){
			addRefFiles(b);
		}else if(arg.equalsIgnoreCase("nt") || arg.equalsIgnoreCase("RefSeq") || arg.equalsIgnoreCase("nr")
				|| arg.equalsIgnoreCase("img") || arg.equalsIgnoreCase("silva") || arg.equalsIgnoreCase("ribo")){
			addRefFiles(arg);
		}else if(a.equals("threads") || a.equals("sketchthreads") || a.equals("t")){
			threads=Integer.parseInt(b);
		}
		
		else if(a.equals("index") || a.equals("makeindex")){
			if(b!=null && "auto".equalsIgnoreCase(b)){
				autoIndex=true;
				makeIndex=true;
			}else{
				autoIndex=false;
				makeIndex=Tools.parseBoolean(b);
			}
		}else if(a.equals("indexsize") || a.equals("indexlimit")){
			SketchIndex.indexLimit=Integer.parseInt(b);
		}
		
		else if(b==null && arg.indexOf('=')<0 && addFileIfNotFound){//if(new File(arg).exists())
			addRefFiles(arg);
		}else{
			return false;
		}
		return true;
	}

	public boolean compare(ArrayList<Sketch> sketches, StringBuilder sb, DisplayParams params, int maxThreads){
		assert(params.postParsed);
		ConcurrentHashMap<Integer, Comparison> map=new ConcurrentHashMap<Integer, Comparison>();
		
		SketchResults[] alca=new SketchResults[sketches.size()];

		if(verbose2){System.err.println("At compare.");}
		
		boolean success=true;
		final CompareBuffer buffer=new CompareBuffer(false);
		AtomicInteger fakeID=new AtomicInteger(minFakeID);
		for(int i=0; i<sketches.size(); i++){
			fakeID.set(minFakeID);
			Sketch a=sketches.get(i);
			
			SketchResults results=processSketch(a, buffer, fakeID, map, params, maxThreads);
			alca[i]=results;
//			System.out.println(a.present);
		}
		
		if(verbose2){System.err.println("Made results.");}
		
		for(int i=0; i<alca.length; i++){
//			Sketch s=sketches.get(i);
			SketchResults results=alca[i];
			sb.append(results.toText(params));
		}
		return success;
	}
	
	private class CompareThread extends Thread {
		
		CompareThread(Sketch a_, ArrayList<Sketch> localRefSketches_, int pid_, int incr_, 
				AtomicInteger fakeID_, ConcurrentHashMap<Integer, Comparison> map_, DisplayParams params_){
			a=a_;
			pid=pid_;
			incr=incr_;
			fakeID=fakeID_;
			map=map_;
			params=params_;
			localRefSketches=localRefSketches_;
			buffer=new CompareBuffer(params.needContamCounts());
			if(buffer.cbs!=null){buffer.cbs.setCapacity(a.length(), 0);}
		}
		
		public void run(){
			assert(a.compareBitSet()==null || buffer.cbs!=null) : (a.compareBitSet()==null)+", "+(buffer.cbs==null); //Unsafe to use a.cbs multithreaded unless atomic
			final AbstractBitSet cbs=(buffer.cbs==null ? a.compareBitSet() : buffer.cbs);
			for(int i=pid; i<localRefSketches.size(); i+=incr){
				Sketch b=localRefSketches.get(i);
				processPair(a, b, buffer, cbs, fakeID, map, params);
			}
		}
		
		final AtomicInteger fakeID;
		final ConcurrentHashMap<Integer, Comparison> map;
		final CompareBuffer buffer;
		final int incr;
		final int pid;
		final Sketch a;
		final DisplayParams params;
		final ArrayList<Sketch> localRefSketches;
		
	}
	
	SketchResults processSketch(Sketch a, CompareBuffer buffer, AtomicInteger fakeID, ConcurrentHashMap<Integer, Comparison> map, DisplayParams params, int maxThreads){
		//		Timer t=new Timer();
		//		t.start("Began query.");
		assert(a.compareBitSet()==null);
		assert(a.indexBitSet()==null);
		
		if(verbose2){System.err.println("At processSketch 1");} //123
		
		a.makeBitSets(params.needContamCounts(), index!=null);
		
		final SketchResults sr;
		if(index!=null){
			sr=index.getSketches(a, params.minHits, params.needContamCounts(), params.contamLevel());
		}else{
			sr=new SketchResults(a, refSketches, null); 
		}
		if(verbose2){System.err.println("At processSketch 2");} //123
		
		if(sr==null || sr.sketchList==null || sr.sketchList.isEmpty()){
			if(verbose2){System.err.println("At processSketch 2.0");} //123
			return sr;
		}
		
		if(verbose2){System.err.println("At processSketch 2.1");} //123
		
		//		t.stop("Got "+hits.size()+" hits.");
		//		t.start();
		//		System.err.println("hits: "+hits.size());
		
		comparisons.getAndAdd(sr.sketchList.size());
		
		if(verbose2){System.err.println("At processSketch 2.2");} //123
		
		if(maxThreads>1 && Shared.threads()>1 && sr.sketchList.size()>31){
			if(verbose2){System.err.println("At processSketch 2.3");} //123
			assert((buffer.cbs==null)==(params.needContamCounts()));
			spawnThreads(a, sr.sketchList, fakeID, map, params, maxThreads);
			if(verbose2){System.err.println("At processSketch 2.4");} //123
		}else{
			if(verbose2){System.err.println("At processSketch 2.5");} //123
			assert(buffer.cbs==null);
			for(Sketch b : sr.sketchList){
				processPair(a, b, buffer, a.compareBitSet(), /*sr.taxHits,*/ fakeID, map, params);
			}
			if(verbose2){System.err.println("At processSketch 2.6");} //123
		}
		if(verbose2){System.err.println("At processSketch 3");} //123
		
		sr.addMap(map, params, buffer);
		
		fakeID.set(minFakeID);
		map.clear();
		if(verbose2){System.err.println("At processSketch 4");} //123
		
		return sr;
	}

	private void spawnThreads(Sketch a, ArrayList<Sketch> refs, AtomicInteger fakeID, 
			ConcurrentHashMap<Integer, Comparison> map, DisplayParams params, int maxThreads){
		final int toSpawn=Tools.max(1, Tools.min((refs.size()+7)/8, threads, maxThreads, Shared.threads()));
		ArrayList<CompareThread> alct=new ArrayList<CompareThread>(toSpawn);
		if(verbose2){System.err.println("At spawnThreads");} //123
		for(int t=0; t<toSpawn; t++){
			alct.add(new CompareThread(a, refs, t, toSpawn, fakeID, map, params));
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
		}
		if(params.needContamCounts()){
			for(CompareThread ct : alct){
				if(ct.buffer.cbs==null){
					assert(AUTOSIZE && index!=null);
					break;
				}
				a.addToBitSet(ct.buffer.cbs);
			}
		}
		alct=null;
	}
	
//	private void writeResults(ArrayList<Comparison> al, Sketch s, StringBuilder sb){
//		sb.append("\nResults for "+s.name()+":\n\n");
//		
//		ArrayList<TaxNode> tnl=new ArrayList<TaxNode>();
//		for(Comparison c : al){
//			formatComparison(c, format, sb, printTax);
//		}
//	}
	
	boolean processPair(Sketch a, Sketch b, CompareBuffer buffer, AbstractBitSet abs, 
			AtomicInteger fakeID, ConcurrentHashMap<Integer, Comparison> map, DisplayParams params){
//		System.err.println("Comparing "+a.name()+" and "+b.name());
		Comparison c=compareOneToOne(a, b, buffer, abs, /*taxHits, params.contamLevel(),*/ params.minHits, params.minWKID, params.minANI, null);
		if(c==null){return false;}
		if(c.taxID()<1){c.taxID=fakeID.getAndIncrement();}
		
//		System.err.println("TID: "+c.taxID()+", "+fakeID);
		
		TaxNode tn=(taxtree==null ? null : taxtree.getNode(b.taxID));
		if(tn!=null){
			c.taxName=tn.name;
			if(tn.level<params.taxLevel){
				tn=taxtree.getNodeAtLevel(b.taxID, params.taxLevel);
			}
		}
		Integer key=(tn==null ? c.taxID : tn.id);
		
		Comparison old=map.get(key);
//		System.err.println("A. Old: "+(old==null ? 0 : old.hits)+", new: "+c.hits);
		if(old!=null && params.compare(old, c)>0){return false;}
		
		old=map.put(key, c);
		while(old!=null && params.compare(old, c)>0){
//			System.err.println("B. Old: "+(old==null ? 0 : old.hits)+", new: "+c.hits);
			c=old;
			old=map.put(key, c);
		}
		return true;
	}
	
	//TODO:  Interestingly, the heap never seems to be created by anything...  not sure what it's for.
	private static Comparison compareOneToOne(final Sketch a, final Sketch b, CompareBuffer buffer, AbstractBitSet abs,
			int minHits, float minWKID, float minANI, Heap<Comparison> heap){
//		assert(heap!=null); //Optional, for testing.
		if(a==b && !compareSelf){return null;}
		final int matches=a.countMatches(b, buffer, abs, true/*!makeIndex || !AUTOSIZE*/ /*, taxHits, contamLevel*/, null, -1);
		assert(matches==buffer.matches);
		if(matches<minHits){return null;}
		
		{
			final int div=buffer.minDivisor();
			final float wkid=matches/(float)div;
			if(wkid<minWKID){return null;}
			//TODO (?)  This is only necessary because of the order of setting minwkid and minani.
			//minWKID can be deterministically determined from minANI so if it is set correctly this can be skipped.
			if(minANI>0){
				final float ani=Comparison.wkidToAni(wkid);
				if(ani<minANI){return null;}
			}
		}
		
		if(heap!=null && !heap.hasRoom() && heap.peek().hits>matches){return null;} //TODO:  Should be based on score
		
//		System.err.print("*");
		Comparison c=new Comparison(buffer, a, b);
		if(heap==null || heap.add(c)){return c;}
		return null;
	}
	
	public void addRefFiles(String a){
		if(a.equalsIgnoreCase("nr")){
			addRefFiles(NR_PATH);
			if(blacklist!=null){blacklist=Blacklist.nrBlacklist();}
			if(defaultParams.dbName==null){defaultParams.dbName="nr";}
			if(!setK){k=9; k2=5;}
		}else if(a.equalsIgnoreCase("nt")){
			addRefFiles(NT_PATH);
			if(blacklist!=null){blacklist=Blacklist.ntBlacklist();}
			if(defaultParams.dbName==null){defaultParams.dbName="nt";}
			if(!setK){k=defaultK; k2=24;}
		}else if(a.equalsIgnoreCase("refseq")){
			addRefFiles(REFSEQ_PATH);
			if(blacklist!=null){blacklist=Blacklist.refseqBlacklist();}
			if(defaultParams.dbName==null){defaultParams.dbName="RefSeq";}
			if(!setK){k=defaultK; k2=24;}
		}else if(a.equalsIgnoreCase("silva")){
//			TaxTree.SILVA_MODE=Tools.parseBoolean(b);
			addRefFiles(SILVA_PATH);
			if(blacklist!=null){blacklist=Blacklist.silvaBlacklist();}
			if(defaultParams.dbName==null){defaultParams.dbName="Silva";}
			if(!setK){k=defaultK; k2=0;}
		}else if(a.equalsIgnoreCase("img")){
			addRefFiles(IMG_PATH);
			if(blacklist!=null){blacklist=Blacklist.imgBlacklist();}
			if(defaultParams.dbName==null){defaultParams.dbName="IMG";}
			if(!setK){k=defaultK; k2=24;}
		}else{
			addFiles(a, refFiles);
		}
	}
	
	private void addFiles(String a, Set<String> list){
		if(a==null){return;}
		File f=new File(a);
		assert(!list.contains(a)) : "Duplicate file "+a;
		
		if(f.exists()){
			list.add(a);
		}else if(a.indexOf(',')>0){
			for(String s : a.split(",")){addFiles(s, list);}
		}else if(a.indexOf('#')>=0 && new File(a.replaceFirst("#", "0")).exists()){
			for(int i=0; true; i++){
				String temp=a.replaceFirst("#", ""+i);
				if(!new File(temp).exists()){break;}
				list.add(temp);
			}
		}else{
			list.add(a);
		}
	}
	
	public void makeIndex(){
		assert(index==null);
		index=new SketchIndex(refSketches);
		index.load();
	}
	
	public void loadReferences(int minKeyOccuranceCount, float minEntropy) {
		makeTool(minKeyOccuranceCount, false);
		refSketches=tool.loadSketches_MT(SketchObject.PER_TAXA, 1f, -1, minEntropy, refFiles);
//		System.err.println("Sketches: "+refSketches.get(0).name());
		if(makeIndex){
			makeIndex();
		}
	}
	
	public void makeTool(int minKeyOccuranceCount, boolean trackCounts){
		if(tool==null){
			tool=new SketchTool(targetSketchSize, minKeyOccuranceCount, trackCounts);
		}
	}
	
	public ArrayList<Sketch> loadSketchesFromString(String sketchString){
		return tool.loadSketchesFromString(sketchString);
	}
	
	public int refFileCount(){return refFiles.size();}
	
	/*--------------------------------------------------------------*/
	
	public SketchIndex index=null;
	public boolean autoIndex=true;
	
	public SketchTool tool=null;
	public ArrayList<Sketch> refSketches=new ArrayList<Sketch>();
	LinkedHashSet<String> refFiles=new LinkedHashSet<String>();
	public int threads=Shared.threads();
	boolean verbose;
	boolean errorState=false;
	AtomicLong comparisons=new AtomicLong(0);
	
}
