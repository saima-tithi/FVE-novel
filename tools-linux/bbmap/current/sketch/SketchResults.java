package sketch;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Map.Entry;
import java.util.concurrent.ConcurrentHashMap;

import fileIO.TextStreamWriter;
import shared.Shared;
import shared.Tools;
import structures.IntHashMap;

public class SketchResults extends SketchObject {
	
	SketchResults(Sketch s){
		sketch=s;
	}
	
	SketchResults(Sketch s, ArrayList<Sketch> sketchList_, int[][] taxHits_){
		sketch=s;
		sketchList=sketchList_;
		taxHits=taxHits_;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	public void addMap(ConcurrentHashMap<Integer, Comparison> map, DisplayParams params, CompareBuffer buffer) {

		if(map.isEmpty()){return;}
		list=addToList(map, params, list);
		
		if((true || params.needContamCounts())){
			recompare(buffer, params);
		}
		if(params.taxFilter!=null){
			int removed=0;
			for(int i=0; i<list.size(); i++){
				Comparison c=list.get(i);
				if((c.taxID()>0 && !params.taxFilter.passesFilter(c.taxID())) || (c.name()!=null && !params.taxFilter.passesFilterByNameOnly(c.name()))){
					list.set(i, null);
					removed++;
				}
			}
			if(removed>0){
				Tools.condenseStrict(list);
			}
		}
	}
	
	public void recompare(CompareBuffer buffer, DisplayParams params){
//		assert(makeIndex || !AUTOSIZE);
		
		assert(!sketch.merged());
		sketch.mergeBitSets();
		
//		System.err.println(sketch.compareBitSet());
//		assert(false) : sketch.compareBitSet().getClass();
		
		for(Comparison c : list){
			c.recompare(buffer, taxHits, params.contamLevel());
		}
		Collections.sort(list, params.comparator);
		Collections.reverse(list);
	}
	
	private static ArrayList<Comparison> addToList(ConcurrentHashMap<Integer, Comparison> map, DisplayParams params, ArrayList<Comparison> old){

//		System.err.println(map.size());
//		System.err.println(map.keySet());
		
		ArrayList<Comparison> al=(old==null ? new ArrayList<Comparison>(map.size()) : old);
		for(Entry<Integer, Comparison> e : map.entrySet()){
			al.add(e.getValue());
		}
		Shared.sort(al, params.comparator);
		Collections.reverse(al);
		while(al.size()>params.maxRecords*2+10){
			al.remove(al.size()-1);
		}
		return al;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Tax Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	public int primaryTax(int level){
		//I have no idea how to implement this...
		IntHashMap map=new IntHashMap();
		assert(false);
		return -1;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Print Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private static String recordBreak="\n"; //"\n\n"
	
	void writeResults(DisplayParams params, TextStreamWriter tsw){
		StringBuilder sb=toText(params);
		tsw.print(sb);
	}

	StringBuilder toText(DisplayParams params){
		assert(params.postParsed);
		final StringBuilder sb=params.queryHeader(sketch);
		if(params.format==3){
			if(list==null || list.isEmpty()){return sb;}
			int idx=0;
			int prevTaxID=0;
			for(Comparison c : list){
				params.formatComparison(c, sb, prevTaxID, sketch);
				prevTaxID=c.taxID();
				idx++;
				if(idx>=params.maxRecords){break;}
			}
		}else{
			sb.append(recordBreak);

			if(list==null || list.isEmpty()){
				sb.append("No hits.\n");
			}else{
				if(params.format==2){sb.append(params.header()).append('\n');}
				int idx=0;
				int prevTaxID=0;
				for(Comparison c : list){
					params.formatComparison(c, sb, prevTaxID, sketch);
					prevTaxID=c.taxID();
					idx++;
					if(idx>=params.maxRecords){break;}
				}
			}
		}
		return sb;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	public final Sketch sketch;
	public ArrayList<Sketch> sketchList;
	public int[][] taxHits;
	public ArrayList<Comparison> list;
	public int totalRecords=0;
	
}
