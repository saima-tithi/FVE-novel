package clump;

import java.io.Serializable;

import hiseq.FlowcellCoordinate;
import shared.Tools;
import stream.KillSwitch;
import stream.Read;
import structures.IntList;

class ReadKey implements Serializable, Comparable<ReadKey> {
	
//	public static ReadKey makeKeyIfNeeded(Read r){
//		if(r.obj==null){
//			return makeKey(r, true);
//		}
//		return (ReadKey)r.obj;
//	}
	
	public static ReadKey makeKey(Read r, boolean setObject){
		assert(r.obj==null);
		try {
			ReadKey rk=new ReadKey(r);
			if(setObject){r.obj=rk;}
			return rk;
		} catch (OutOfMemoryError e) {
			KillSwitch.memKill(e);
			throw new RuntimeException(e);
		}
	}
	
	private ReadKey(Read r){
		this(r, 0, 0, true);
	}
	
	private ReadKey(Read r, long kmer_, int position_, boolean plus_){
		kmer=kmer_;
		position=position_;
		clump=null;
		assert(!r.swapped());
		flipped=false;
		kmerMinusStrand=false;
		
		if(Clump.opticalOnly){
			FlowcellCoordinate fc=FlowcellCoordinate.getFC();
			fc.setFrom(r.id);
			lane=fc.lane;
			tile=fc.tile;
			x=fc.x;
			y=fc.y;
		}
		
	}
	
	protected ReadKey(){}
	
	public void set(long kmer_, int position_, boolean minus_){
		setKmer(kmer_);
		setPosition(position_);
//		setClump(null);
		kmerMinusStrand=minus_;
	}
	
	private long setKmer(long kmer_){
		return kmer=kmer_;
	}
	
	private int setPosition(int position_){
		return position=position_;
	}
	
	public Clump setClump(Clump clump_){
		return clump=clump_;
	}
	
	private boolean setFlipped(boolean flipped_){
		assert(flipped!=flipped_);
		return flipped=flipped_;
	}
	
	public void clear(){
		setKmer(0);
		setPosition(0);
		setClump(null);
		kmerMinusStrand=false;
	}
	
	public void flip(Read r, int k){
		assert(r.swapped()==flipped);
		r.reverseComplement();
		r.setSwapped(!r.swapped());
		setFlipped(!flipped);
		if(r.length()>=k){setPosition(r.length()-position+k-2);}
		assert(r.swapped()==flipped);
	}
	
	@Override
	public int compareTo(ReadKey b){
		if(kmer!=b.kmer){return kmer>b.kmer ? -1 : 1;} //Bigger kmers first...
		if(kmerMinusStrand!=b.kmerMinusStrand){return kmerMinusStrand ? 1 : -1;}
		if(position!=b.position){return position<b.position ? 1 : -1;}
		if(Clump.opticalOnly){
			if(lane!=b.lane){return lane-b.lane;}
			if(!spanTiles()){
				if(tile!=b.tile){return tile-b.tile;}
			}
			if(Clump.sortYEarly()){ //Improves speed slightly
				if(y!=b.y){return y-b.y;}
			}
		}
		return 0;
	}
	
	@Override
	public boolean equals(Object b){
		return equals((ReadKey)b);
	}
	
	/** Note that this is different than compareTo()==0
	 * That's to prevent sortYEarly comparison making things unequal.
	 * @param b
	 * @return
	 */
	public boolean equals(ReadKey b){
		if(b==null){return false;}
		if(kmer!=b.kmer || kmerMinusStrand!=b.kmerMinusStrand || position!=b.position){return false;}
		if(Clump.opticalOnly){
			if(lane!=b.lane){return false;}
			if(!spanTiles()){
				if(tile!=b.tile){return false;}
			}
		}
		return true;
	}
	
	public String toString(){
		return position+","+(kmerMinusStrand ? ",t" : ",f")+","+kmer+"\t"+lane+","+tile+","+x+","+y;
	}
	
	public float distance(ReadKey rkb){
		if(lane!=rkb.lane){return FlowcellCoordinate.big;}
		
		long a=Tools.absdif(x, rkb.x), b=Tools.absdif(y, rkb.y);
		if(tile!=rkb.tile){
			if(spanTiles()){
				if(spanAdjacentOnly && Tools.absdif(tile, rkb.tile)>1){return FlowcellCoordinate.big;}
				if(spanTilesX && spanTilesY){
					return Tools.min(a, b);
				}else if(spanTilesX){
					return a;
				}else{
					return b;
				}
			}else{
				return FlowcellCoordinate.big;
			}
		}
		return (float)Math.sqrt(a*a+b*b);
	}
	
	public boolean near(ReadKey rkb, float dist){
		return distance(rkb)<dist;
	}
	
	public boolean nearXY(ReadKey rkb, float dist){
		if(lane!=rkb.lane){return false;}
		
		long a=Tools.absdif(x, rkb.x), b=Tools.absdif(y, rkb.y);
		return Tools.min(a,b)<=dist;
	}
	
	public int flag;
	
	public long kmer;
	public int position;
	public boolean flipped;
	public boolean kmerMinusStrand;
	public Clump clump;
	public IntList vars;
	
	public int lane, tile, x, y;
	
	public static final int overhead=overhead();
	private static int overhead(){
		return 16+ //self
				1*(8)+ //kmer
				1*(4)+ //int fields
				2*(4)+ //booleans
				2*(8)+ //pointers
				4*(4); //flowcell coordinate
	}
	
	public static boolean spanTiles(){return spanTilesX || spanTilesY;}
	public static boolean spanTilesX=false;
	public static boolean spanTilesY=false;
	public static boolean spanAdjacentOnly=false;
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	
}
