package tax;

import java.io.Serializable;
import java.util.Comparator;

import shared.Tools;

/**
 * Represents a taxonomic identifier, such as a specific genus.
 * Includes the name, NCBI numeric id, parent id, and taxonomic level.
 * @author Brian Bushnell
 * @date Mar 6, 2015
 *
 */
public class TaxNode implements Serializable{
	
	private static final long serialVersionUID = 315375359526248450L;

	public TaxNode(int id_, String name_){
		this(id_, -1, -1, -1, name_);
	}
	
	public TaxNode(int id_, int parent_, int level_, int levelExtended_, String name_){
		id=id_;
		pid=parent_;
		level=level_;
		levelExtended=levelExtended_;
		setOriginalLevel(levelExtended);
		name=name_;
	}

	/**
	 * @param split
	 * @param i
	 * @return
	 */
	public boolean matchesName(String[] split, int idx, TaxTree tree) {
		if(idx<0){return true;}
		if(!split[idx].equalsIgnoreCase(name)){return false;}
		return tree.getNode(pid).matchesName(split, idx-1, tree);
	}
	
	public String toString(){
		return "("+id+","+pid+","+countRaw+","+countSum+",'"+levelStringExtended(false)+"',"+(canonical() ? "T" : "F")+",'"+name+"')";
	}
	
	public boolean equals(TaxNode b){
		if(id!=b.id || pid!=b.pid || levelExtended!=b.levelExtended || flag!=b.flag){return false;}
		if(name==b.name){return true;}
		if((name==null) != (b.name==null)){return false;}
		return name.equals(b.name);
	}
	
	public long incrementRaw(long amt){
		if(amt==0){return countRaw;}
		if(verbose){System.err.println("incrementRaw("+amt+") node: "+this);}
		countRaw+=amt;
		assert(countRaw>=0) : "Overflow! "+countRaw+", "+amt;
		return countRaw;
	}
	
	public long incrementSum(long amt){
		if(amt==0){return countSum;}
		if(verbose){System.err.println("incrementSum("+amt+") node: "+this);}
		countSum+=amt;
		assert(countSum>=0) : "Overflow! "+countSum+", "+amt;
		return countSum;
	}
	
	public boolean isSimple(){
		return levelExtended!=TaxTree.NO_RANK_E && levelExtended==TaxTree.levelToExtended(level);
	}
	
//	public String levelString(){return level<0 ? "unknown" : TaxTree.levelToString(level);}
	
	public String levelStringExtended(boolean original){
		int x=(original ? originalLevel() : levelExtended);
		return x<0 ? "unknown" : TaxTree.levelToStringExtended(x);
	}

	public String levelToStringShort() {return level<0 ? "x" : TaxTree.levelToStringShort(level);}
	
	public static class CountComparator implements Comparator<TaxNode>{
		
		@Override
		public int compare(TaxNode a, TaxNode b) {
			long x=b.countSum-a.countSum;
//			System.err.println("x="+x+" -> "+Tools.longToInt(x));
			if(x!=0){return Tools.longToInt(x);}
			return a.levelExtended==b.levelExtended ? a.id-b.id : a.levelExtended-b.levelExtended;
		}
		
	}
	
	public void setCanonical(boolean b){
		if(b){flag=flag|CANON_MASK;}
		else{flag=flag&~CANON_MASK;}
	}
	
	public void setOriginalLevel(int x){
		flag=(flag&~ORIGINAL_LEVEL_MASK)|(x&ORIGINAL_LEVEL_MASK);
	}
	
	public boolean canonical(){
		return (flag&CANON_MASK)==CANON_MASK;
	}
	
	public boolean levelChanged(){
		return originalLevel()!=levelExtended;
	}
	
	public int originalLevel(){
		int x=(int)(flag&ORIGINAL_LEVEL_MASK);
		return x==ORIGINAL_LEVEL_MASK ? -1 : x;
	}
	
	public boolean cellularOrganisms(){
		return id==TaxTree.CELLULAR_ORGANISMS;
	}
	
	public static final long ORIGINAL_LEVEL_MASK=63;
	public static final long CANON_MASK=64;
	
	@Override
	public final int hashCode(){return id;}

	public final int id;
	public final String name;
	public int pid;
	public int level;
	public int levelExtended;
	
	private long flag=0;

	public long countRaw=0;
	public long countSum=0;
	
	public static final boolean verbose=false;
	public static final CountComparator countComparator=new CountComparator();
	
	
}
