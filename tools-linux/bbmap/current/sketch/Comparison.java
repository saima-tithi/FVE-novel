package sketch;

import java.util.Arrays;
import java.util.Comparator;
import java.util.Locale;

import shared.Tools;

public final class Comparison extends SketchObject implements Comparable<Comparison> {
	
	/*--------------------------------------------------------------*/
	/*----------------         Constructors         ----------------*/
	/*--------------------------------------------------------------*/
	
	public Comparison(){}
	
	public Comparison(CompareBuffer buffer){
		this(buffer, null, null);
	}
	
	public Comparison(Sketch a_, Sketch b_){
		this(null, a_, b_);
	}
	
	public Comparison(CompareBuffer buffer, Sketch a_, Sketch b_){
		
		a=a_;
		b=b_;
		
		if(buffer!=null){setFrom(buffer);}
		
		if(b!=null){
			taxName=b.taxName();
			taxID=b.taxID;
		}

//		System.err.println(this);
//		System.err.println(b.present);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Mutators           ----------------*/
	/*--------------------------------------------------------------*/
	
	public void setFrom(CompareBuffer buffer){
		hits=buffer.matches;
		multiHits=buffer.multiMatches;
		unique2=buffer.unique2;
		unique3=buffer.unique3;
		noHits=buffer.noHits;

		contamHits=buffer.contamHits;
		contam2Hits=buffer.contam2Hits;
		multiContamHits=buffer.multiContamHits;
		
		refDivisor=buffer.refDivisor;
		queryDivisor=buffer.queryDivisor;
		
		refSize=buffer.refSize;
		querySize=buffer.querySize;

		depth=buffer.depth();
		depth2=buffer.depth2();
//		volume=volume0();
		score=score0();
	}
	
	public void recompare(CompareBuffer buffer, int[][] taxHits, int contamLevel){
		
//		for(int[] row : taxHits){
//			if(row!=null){
//				System.err.println(Arrays.toString(row));
//			}
//		}//Tested; correctly indicates most rows have octopus but some have other things.
		
		assert(a.merged());
//		int oldContam2=contam2Hits;
		int x=a.countMatches(b, buffer, a.compareBitSet(), false, taxHits, contamLevel);
		assert(x==hits);
		setFrom(buffer);
//		contam2Hits=oldContam2;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	public boolean equals(Object b){
		if(b==null || b.getClass()!=this.getClass()){return false;}
		return scoreComparator.compare(this, (Comparison)b)==0;
	}
	
	//WKID
	public float idMinDivisor(){
		return hits/(float)minDivisor();
	}
	
	//KID
	public float idMaxDivisor(){
		return hits/(float)maxDivisor();
	}
	
	public float idQueryDivisor(){
		return hits/(float)(Tools.max(1, refDivisor));
	}
	
	public float idRefDivisor(){
		return hits/(float)(Tools.max(1, refDivisor));
	}
	
	public float completeness(){
		float complt=Tools.min(1, (queryDivisor-contamHits)/(float)Tools.max(1, refDivisor));
		return complt;
//		float c2=hits/(float)refDivisor;
//		assert(queryDivisor-contamHits>=hits);
//		assert(c1>=c2);
//		System.err.println(hits+", "+contamHits+", "+refDivisor+", "+queryDivisor+", "+c1+", "+c2);
//		return Tools.max(c1, c2);
//		float kid=idMaxDivisor(), wkid=idMinDivisor();
//		return kid==0 ? 0 : kid/wkid;
	}
	
	public float contamFraction(){
		return Tools.min(1, contamHits/(float)Tools.max(1, queryDivisor));
	}
	
	public float contam2Fraction(){
		return Tools.min(1, contam2Hits/(float)Tools.max(1, queryDivisor));
	}
	
	public float uContamFraction() {
		int uContamHits=contamHits-multiContamHits;
		return Tools.min(1, uContamHits/(float)Tools.max(1, queryDivisor));
	}
	
	public float ani(){
		if(hits<1){return 0;}
		
		double wkid=idMinDivisor();
		return wkidToAni(wkid);

//		final float rID=hits/(float)(refDivisor);
//		final float qID=hits/(float)(queryDivisor-contamHits);
//		final float wkid2=Tools.max(qID, rID);
//		final float ani=wkidToAni(wkid2);
//		
////		System.err.println("rid: "+wkidToAni(rID)+", qid: "+wkidToAni(qID)+", qid2: "+wkidToAni(hits/(float)(queryDivisor)));
//		
//		return ani;
	}
	
	int minDivisor(){return Tools.max(1, Tools.min(refDivisor, queryDivisor));}
	int maxDivisor(){return Tools.max(1, refDivisor, queryDivisor);}
	
	private float score0(){
		long est=useSizeEstimate ? genomeSizeEstimate() : genomeSizeKmers();
		float wkid=idMinDivisor();
		float kid=idMaxDivisor();
		float complt=completeness();
		float contam=contamFraction();
		return (float)(Math.sqrt(80*(40000+hits+uHits())*(wkid*kid*Math.pow(est*complt, 0.2)*(1-contam*0.1)))+0.1);
	}
	
	public String scoreS(){
		float x=score;
		return format3(x);
	}
	
	public String depthS(boolean observedToActual){
		float x=depth;
		if(observedToActual){x=(float)Tools.observedToActualCoverage(x);}
		return format3(x);
	}
	
	public String depth2S(boolean observedToActual){
		float x=depth2;
		if(observedToActual){
			x=(float)(Tools.observedToActualCoverage(depth)*(depth2/depth));
		}
		return format3(x);
	}
	
	public String volumeS(){
		double x=volume()*0.001;
		return format3(x);
	}
	
	String format3(double x){
		if(x>=999.95){
			return(""+(long)Math.round(x));
		}else if(x>=99.995){
			return String.format(Locale.ROOT, "%.1f", x);
		}else if(x>=9.9995){
			return String.format(Locale.ROOT, "%.2f", x);
		}
		return String.format(Locale.ROOT, "%.3f", x);
	}
	
	private float volume(){
		return Tools.max(1f, depth)*hits;
	}
	
	public String toString(){
		return "hits="+hits+", refDivisor="+refDivisor+", queryDivisor="+queryDivisor+", refSize="+refSize+", querySize="+querySize+
				", contamHits="+contamHits+", contam2Hits="+contam2Hits+", multiContamHits="+multiContamHits+", depth="+depth+", depth2="+depth2+", volume="+volume()+
				", hits="+hits+", multiHits="+multiHits+", unique2="+unique2+", unique3="+unique3+", noHits="+noHits+", taxID="+taxID+", taxName="+taxName;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Getters           ----------------*/
	/*--------------------------------------------------------------*/

	public String name(){return taxName!=null ? taxName : name0()!=null ? name0() : fname();}
	public String taxName(){return taxName;}
	String name0(){return b.name0();}
	String fname(){return b.fname();}

//	public int taxID(){return b.taxID<minFakeID ? b.taxID : 0;}
	public int taxID(){return (taxID<minFakeID && taxID>=0) ? taxID : 0;}
	long imgID(){return (b.imgID>0 ? b.imgID : 0);}
	
	long genomeSizeBases(){return b.genomeSizeBases;}
	long genomeSizeKmers(){return b.genomeSizeKmers;}
	long genomeSequences(){return b.genomeSequences;}
	long genomeSizeEstimate(){return b.genomeSizeEstimate();}

	public int uHits() {return hits-multiHits;}
	
	/*--------------------------------------------------------------*/
	/*----------------         Comparators          ----------------*/
	/*--------------------------------------------------------------*/
	
	
	
	static class ScoreComparator implements Comparator<Comparison>{

		@Override
		public int compare(Comparison a, Comparison b) {
			{
				float pa=a.score, pb=b.score;
				if(pa>pb){
					return 1;
				}else if (pa<pb){
					return -1;
				}
			}
			
			int x=a.hits-b.hits;
			if(x!=0){return x;}
			x=b.minDivisor()-a.minDivisor();
			if(x!=0){return x;}
			x=b.maxDivisor()-a.maxDivisor();
			if(x!=0){return x;}
			x=b.refDivisor-a.refDivisor;
			if(x!=0){return x;}
			x=a.taxID()-b.taxID();
			if(x!=0){return x;}
			if(a.name0()!=null && b.name0()!=null){
				return a.name0().compareTo(b.name0());
			}
			if(a.taxName()!=null && b.taxName()!=null){
				return a.taxName().compareTo(b.taxName());
			}
			return 0;
		}
		
		public String toString(){return "sortByScore";}
		
	}
	
	static class DepthComparator implements Comparator<Comparison>{

		@Override
		public int compare(Comparison a, Comparison b) {
			final float da=Tools.max(0.1f, a.depth-0.5f), db=Tools.max(0.1f, b.depth-0.5f);
			final float sa, sb;
			if(sqrt){
				sa=da*(float)Math.sqrt(a.score);
				sb=db*(float)Math.sqrt(b.score);
			}else{
				sa=da*a.score;
				sb=db*b.score;
			}
			return sa>sb ? 1 : sa<sb ? -1 : scoreComparator.compare(a, b);
		}
		
		public String toString(){return "sortByDepth";}
		
	}
	
	static class Depth2Comparator implements Comparator<Comparison>{

		@Override
		public int compare(Comparison a, Comparison b) {
			final float da=Tools.max(0.1f, a.depth2-0.8f), db=Tools.max(0.1f, b.depth2-0.8f);
			final float sa, sb;
			if(sqrt){
				sa=da*(float)Math.sqrt(a.score);
				sb=db*(float)Math.sqrt(b.score);
			}else{
				sa=da*a.score;
				sb=db*b.score;
			}
			return sa>sb ? 1 : sa<sb ? -1 : scoreComparator.compare(a, b);
		}
		
		public String toString(){return "sortByDepth2";}
		
	}
	
	static class VolumeComparator implements Comparator<Comparison>{

		@Override
		public int compare(Comparison a, Comparison b) {
			final float da=a.volume(), db=b.volume();
			final float sa, sb;
			if(sqrt){
				sa=da*(float)Math.sqrt(a.score);
				sb=db*(float)Math.sqrt(b.score);
			}else{
				sa=da*a.score;
				sb=db*b.score;
			}
			return sa>sb ? 1 : sa<sb ? -1 : scoreComparator.compare(a, b);
		}
		
		public String toString(){return "sortByVolume";}
		
	}
	
	@Override
	public int compareTo(Comparison b) {
		assert(false) : "Please use comparators instead.";
		return scoreComparator.compare(this, b);
	}

	public static final ScoreComparator scoreComparator=new ScoreComparator();
	public static final DepthComparator depthComparator=new DepthComparator();
	public static final Depth2Comparator depth2Comparator=new Depth2Comparator();
	public static final VolumeComparator volumeComparator=new VolumeComparator();
	private static final boolean sqrt=false;
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	public Sketch a, b;

	String taxName;
	int taxID;
	
	int hits;
	int multiHits;
	int unique2;
	int unique3;
	int noHits;

	float depth;
	float depth2;
	float score;

	int contamHits;
	int contam2Hits;
	int multiContamHits;
	
	int refDivisor;
	int queryDivisor;
	
	int refSize;
	int querySize;
	
}
