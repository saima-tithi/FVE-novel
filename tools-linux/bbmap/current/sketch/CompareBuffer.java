package sketch;

import shared.Tools;
import structures.AbstractBitSet;

public class CompareBuffer {
	
	public CompareBuffer(boolean makeBS){
		if(makeBS){
			cbs=AbstractBitSet.make(0, SketchObject.bitSetBits);
		}else{
			cbs=null;
		}
	}
	
	void set(final int matches_, final int multiMatches_, final int unique2_, final int unique3_, final int noHits_,
			final int contamHits_, final int contam2Hits_, final int multiContamHits_,
			final int queryDivisor_, final int refDivisor_, final int querySize_, final int refSize_, final long depthSum_, final double depthSum2_){
		matches=matches_;
		multiMatches=multiMatches_;
		unique2=unique2_;
		unique3=unique3_;
		noHits=noHits_;
		
		contamHits=contamHits_;
		contam2Hits=contam2Hits_;
		multiContamHits=multiContamHits_;
		
		queryDivisor=queryDivisor_;
		refDivisor=refDivisor_;
		
		querySize=querySize_;
		refSize=refSize_;

		depthSum=depthSum_;
		depthSum2=(float)depthSum2_;
	}
	
	void clear(){
		matches=multiMatches=0;
		unique2=unique3=noHits=0;
		contamHits=contam2Hits=multiContamHits=0;
		refDivisor=queryDivisor=0;
		refSize=querySize=0;
		depthSum=0;
		depthSum2=0;
	}
	
	float depth(){
		return depthSum<1 ? 0 : depthSum/Tools.max(1.0f, matches);
	}
	
	float depth2(){
		return depthSum2<=0 ? 0 : depthSum2/Tools.max(1.0f, matches);
	}

	int minDivisor(){return Tools.max(1, Tools.min(queryDivisor, refDivisor));}
	int maxDivisor(){return Tools.max(1, queryDivisor, refDivisor);}
	int minSize(){return Tools.max(1, Tools.min(querySize, refSize));}
	int maxSize(){return Tools.max(1, querySize, refSize);}

	int uniqueMatches(){return matches-multiMatches;}
	int uniqueContamHits(){return contamHits-multiContamHits;}
	
	int matches;
	int multiMatches;
	int noHits;
	int unique2;
	int unique3;

	int contamHits;
	int contam2Hits;
	int multiContamHits;
	
	int queryDivisor;
	int refDivisor;
	
	int querySize;
	int refSize;

	long depthSum;
	float depthSum2;

	public final AbstractBitSet cbs; //Only for comparisons, not index
	
}
