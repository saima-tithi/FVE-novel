package sketch;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import dna.AminoAcid;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.EntropyTracker;
import tax.TaxTree;

public class SketchObject {
	
	/*--------------------------------------------------------------*/
	/*----------------           Parsing            ----------------*/
	/*--------------------------------------------------------------*/
	
	public static void main(String[] args){
		Timer t=new Timer();
		for(int i=0; i<1000000; i++){
			wkidToAniExact(0.5, 31, 24, 0.000005);
		}
		t.stop();
		System.out.println(t);
	}
	
	public static boolean parseSketchFlags(String arg, String a, String b){
		
		if(parseCoding(arg, a, b)){
			//
		}
		
		else if(a.equalsIgnoreCase("k")){
			if(b.indexOf(',')>=0){
				String[] split=b.split(",");
				assert(split.length==2) : "Bad argument "+arg;
				int x=Integer.parseInt(split[0]);
				int y=Integer.parseInt(split[1]);
				k=Tools.max(x, y);
				k2=Tools.min(x, y);
			}else{
				k=Integer.parseInt(b);
				k2=0;
			}
			setK=true;
		}else if(a.equalsIgnoreCase("k2")){
			k2=Integer.parseInt(b);
			setK=true;
		}else if(a.equalsIgnoreCase("rcomp")){
			rcomp=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("minfakeid")){
			minFakeID=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("hashNames")){
			hashNames=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("skipNonCanonical")){
			skipNonCanonical=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("useSizeEstimateInScore") || a.equalsIgnoreCase("useSizeEstimate")){
			useSizeEstimate=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("minfakeid")){
			minFakeID=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("blacklist")){
			blacklist=b;
			if(b!=null){
				if("nt".equalsIgnoreCase(b) && !new File(b).exists()){
					blacklist=Blacklist.ntBlacklist();
				}else if(("silva".equalsIgnoreCase(b) || "ribo".equalsIgnoreCase(b)) && !new File(b).exists()){
					blacklist=Blacklist.silvaBlacklist();
				}else if("refseq".equalsIgnoreCase(b) && !new File(b).exists()){
					blacklist=Blacklist.refseqBlacklist();
				}else if("img".equalsIgnoreCase(b) && !new File(b).exists()){
					blacklist=Blacklist.imgBlacklist();
				}
			}
		}else if(a.equalsIgnoreCase("whitelist") || a.equalsIgnoreCase("usewhitelist")){
			useWhitelist=Tools.parseBoolean(b);
		}
		
		else if(a.equalsIgnoreCase("sampleseed")){
			sampleseed=Long.parseLong(b);
		}
		
		else if(a.equalsIgnoreCase("size") || a.equalsIgnoreCase("sketchsize")){
			if("auto".equalsIgnoreCase(b)){
				AUTOSIZE=true;
			}else{
				AUTOSIZE=false;
				targetSketchSize=(int)Tools.parseKMG(b);
			}
		}else if(a.equalsIgnoreCase("maxfraction") || a.equalsIgnoreCase("maxgenomefraction") || a.equalsIgnoreCase("mgf")){
			maxGenomeFraction=Float.parseFloat(b);
		}else if(a.equalsIgnoreCase("smallsketchsize")){
			smallSketchSize=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("minSketchSize")){
			minSketchSize=Tools.max(1, Integer.parseInt(b));
		}else if(a.equalsIgnoreCase("autosize")){
			AUTOSIZE=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("autosizefactor") || a.equalsIgnoreCase("autosizefraction") || a.equalsIgnoreCase("autosizemult") || a.equalsIgnoreCase("sizemult")){
			AUTOSIZE_FACTOR=Float.parseFloat(b);
		}else if(a.equalsIgnoreCase("maxGenomeFractionSmall")){
			maxGenomeFractionSmall=Float.parseFloat(b);
		}else if(a.equalsIgnoreCase("keyFraction")){
			double d=Double.parseDouble(b);
			setKeyFraction(d);
		}
		
		else if(a.equalsIgnoreCase("intmap")){
			SketchIndex.useIntMap=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("intmapbinary")){
			boolean x=Tools.parseBoolean(b);
			SketchIndex.useIntMapBinary=x;
			if(x){SketchIndex.useIntMap=true;}
		}else if(a.equalsIgnoreCase("intmapsize")){
			SketchIndex.intMapSize=(int)Tools.parseKMG(b);
		}else if(a.equalsIgnoreCase("bitsetbits")){
			bitSetBits=Integer.parseInt(b);
		}
		
//		else if(a.equalsIgnoreCase("minkmercount") || a.equalsIgnoreCase("minkeycount")){
//			minKeyOccuranceCount=Integer.parseInt(b);
//		}
		else if(a.equalsIgnoreCase("sketchHeapFactor")){
			sketchHeapFactor=Integer.parseInt(b);
		}
		
		else if(a.equalsIgnoreCase("killok")){
			KILL_OK=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("exactani")){
			EXACT_ANI=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("hashSeed")){
			hashSeed=Long.parseLong(b);
			remakeCodes(hashSeed);
		}
		
		else if(a.equalsIgnoreCase("minprob") || a.equalsIgnoreCase("pfilter")){
			minProb=(float)Double.parseDouble(b);
		}else if(a.equalsIgnoreCase("minQual") || a.equalsIgnoreCase("minq")){
			minQual=Byte.parseByte(b);
		}
		
		else if(a.equalsIgnoreCase("forceDisableMultithreadedFastq")){
			forceDisableMultithreadedFastq=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("merge") || a.equalsIgnoreCase("mergePairs") || a.equalsIgnoreCase("ecco")){
			mergePairs=Tools.parseBoolean(b);
		}
		
		else if(a.equalsIgnoreCase("verifyentropy")){
			EntropyTracker.verify=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("entropyK")){
			entropyK=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("fastentropy") || a.equalsIgnoreCase("fentropy")){
			if(Tools.parseBoolean(b)){EntropyTracker.speed=EntropyTracker.FAST;}
		}else if(a.equalsIgnoreCase("mediumentropy") || a.equalsIgnoreCase("mentropy")){
			if(Tools.parseBoolean(b)){EntropyTracker.speed=EntropyTracker.MEDIUM;}
		}else if(a.equalsIgnoreCase("slowentropy") || a.equalsIgnoreCase("sentropy")){
			if(Tools.parseBoolean(b)){EntropyTracker.speed=EntropyTracker.SLOW;}
		}else if(a.equalsIgnoreCase("superslowentropy") || a.equalsIgnoreCase("ssentropy")){
			if(Tools.parseBoolean(b)){EntropyTracker.speed=EntropyTracker.SUPERSLOW;}
		}
		
//		else if(a.equalsIgnoreCase("minHashValue")){
//			minHashValue=Tools.max(1, Long.parseLong(b));
//		}
		
		else{
			return false;
		}
		return true;
	}
	
	private static boolean parseCoding(String arg, String a, String b){
		if(a.equalsIgnoreCase("deltaout") || a.equalsIgnoreCase("delta")){
			deltaOut=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("a33") || a.equalsIgnoreCase("a48")){
			boolean x=Tools.parseBoolean(b);
			if(x){CODING=A48;}
			else if(CODING==A48){CODING=HEX;}
		}else if(a.equalsIgnoreCase("hex")){
			boolean x=Tools.parseBoolean(b);
			if(x){CODING=HEX;}
			else if(CODING==HEX){CODING=A48;}
		}else if(a.equalsIgnoreCase("raw")){
			boolean x=Tools.parseBoolean(b);
			if(x){CODING=RAW;}
			else if(CODING==RAW){CODING=A48;}
		}else{
			return false;
		}
		return true;
	}
	
	static int parseMode(String arg, String a, String b){
		int mode_=-1;
		if(a.equalsIgnoreCase("mode")){
			if(b!=null && Character.isDigit(b.charAt(0))){
				mode_=Integer.parseInt(b);
			}else if(b.equalsIgnoreCase("single") || b.equalsIgnoreCase("onesketch")){
				mode_=ONE_SKETCH;
			}else if(b.equalsIgnoreCase("taxa") || b.equalsIgnoreCase("pertaxa")){
				mode_=PER_TAXA;
			}else if(b.equalsIgnoreCase("sequence") || b.equalsIgnoreCase("persequence")){
				mode_=PER_SEQUENCE;
			}else if(b.equalsIgnoreCase("header") || b.equalsIgnoreCase("perheader")){
				mode_=PER_HEADER;
			}else if(b.equalsIgnoreCase("img")){
				mode_=PER_IMG;
			}else if(Character.isDigit(b.charAt(0))){
				mode_=Integer.parseInt(b);
			}
		}else if(a.equalsIgnoreCase("single") || a.equalsIgnoreCase("onesketch")){
			mode_=ONE_SKETCH;
		}else if(a.equalsIgnoreCase("taxa") || a.equalsIgnoreCase("pertaxa")){
			mode_=PER_TAXA;
		}else if(a.equalsIgnoreCase("sequence") || a.equalsIgnoreCase("persequence")){
			mode_=PER_SEQUENCE;
		}else if(a.equalsIgnoreCase("header") || a.equalsIgnoreCase("perheader")){
			mode_=PER_HEADER;
		}
		return mode_;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	static synchronized void setTaxtree(String taxTreeFile){
		if(taxTreeFile==null){
			taxtree=null;
			return;
		}
		if(treefile!=null){
			assert(!treefile.equals(taxTreeFile));
			if(treefile.equals(taxTreeFile)){return;}
			treefile=taxTreeFile;
		}
		taxtree=TaxTree.loadTaxTree(taxTreeFile, System.err, hashNames);
	}
	
	public static void reset(){
		postparsed=false;
		blacklist=null;
		useWhitelist=false;
		defaultParams=new DisplayParams();
	}
	
	public static void postParse(){
		if(postparsed){return;}
		postparsed=true;
		
		if(defaultParams.amino()){amino=true;}
		if(amino){Shared.AMINO_IN=true;}
		
		if(amino8){
			AminoAcid.AMINO_SHIFT=3;
			System.err.println("Set AMINO_SHIFT to "+AminoAcid.AMINO_SHIFT);
		}
		
		if(amino){
			int oldk=k;
			int oldk2=k2;
			k=Tools.min(k, 63/AminoAcid.AMINO_SHIFT);
			k2=Tools.min(k2, k);
			if(k2==k){k2=0;}
			rcomp=false;
			if(k!=oldk || k2!=oldk2){System.err.println("Set k to "+k+(k2>0 ? ","+k2 : ""));}
//			AUTOSIZE_FACTOR=(AUTOSIZE_FACTOR*defaultK)/k;
//			System.err.println("Set AUTOSIZE_FACTOR to "+AUTOSIZE_FACTOR);
		}
		
		if(amino){
			bitsPerBase=4;
		}else{
			bitsPerBase=2;
		}
		basesPerCycle=bitsPerCycle/bitsPerBase;
		hashCycles=((k*bitsPerBase)+bitsPerCycle-1)/bitsPerCycle;
		
		if(k2>0){
			assert(k2<k) : "k2 must be less than k.";
//			assert(k2%basesPerCycle==0) : "k2 must be a multiple of "+basesPerCycle; //No longer required!  Still recommended for speed though.
			
			k2mask=~((-1L)<<(bitsPerBase*k2));
			k2submask=~((-1L)<<(bitsPerBase*(k2%basesPerCycle)));
			hashCycles2=k2/basesPerCycle;
		}
		codeMax=hashCycles*codeIncrement;
		codeMax2=hashCycles2*codeIncrement;
		
//		hasher=k2<1 ? new Hasher1() : k2submask==0 ? new Hasher2() : new Hasher3();
		
		if(blacklist!=null){Blacklist.addFiles(blacklist);}
		
//		System.err.println("amino="+amino+", k="+k+", k2="+k2+", bitsPerCycle="+bitsPerCycle+", bitsPerBase="+bitsPerBase+
//				", basesPerCycle="+basesPerCycle+", hashCycles="+hashCycles+", k2mask="+k2mask+", k2submask="+k2submask+", hashCycles2="+hashCycles2+
//				", codeMax="+codeMax+", codeMax2="+codeMax2);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Hashing            ----------------*/
	/*--------------------------------------------------------------*/
	
	private static void remakeCodes(long hashSeed){
		long[][] codes2=makeCodes(maxCycles, codeIncrement, hashSeed, false);
		long[] codes1D2=makeCodes1D(maxCycles, codeIncrement, hashSeed, false);
		for(int i=0; i<codes2.length; i++){
			for(int j=0; j<codes2[i].length; j++){
				codes[i][j]=codes2[i][j];
			}
		}
		for(int i=0; i<codes1D2.length; i++){
				codes1D[i]=codes1D2[i];
		}
	}
	
	public static long[][] makeCodes(int symbols, int modes, long seed, boolean positive){
		Random randy=(seed>=0 ? new Random(seed) : new Random());
		long mask=positive ? Long.MAX_VALUE : -1L;
		long[][] r=new long[symbols][modes];
		for(int i=0; i<symbols; i++){
			for(int j=0; j<modes; j++){
				r[i][j]=randy.nextLong()&mask;
			}
		}
		return r;
	}
	
	public static long[] makeCodes1D(int symbols, int modes, long seed, boolean positive){
		Random randy=(seed>=0 ? new Random(seed) : new Random());
		long mask=positive ? Long.MAX_VALUE : -1L;
		long[] r=new long[symbols*modes];
		for(int i=0; i<r.length; i++){
			r[i]=randy.nextLong()&mask;
		}
		return r;
	}
	
	public static final long hash(long kmer){
		assert(postparsed);
//		return hasher.hash_inner(kmer); //Slow
//		return k2<1 ? hash1(kmer) : hash2(kmer);
		return k2<1 ? hash1(kmer) : k2submask==0 ? hash2(kmer) : hash3(kmer);
	}
	
//	/**
//	 * Generates a hash code from a kmer.
//	 * @param kmer Kmer to hash
//	 * @return Hash code
//	 */
//	private static final long hash1(long kmer){
//		long code=kmer;
//		for(int i=0; i<hashCycles; i++){
//			int x=(int)(kmer&cycleMask);
//			kmer>>=bitsPerCycle;
//			code^=codes[i][x];
////			code^=Long.rotateRight(codes[i][x], (int)(code&3)); //Might be better.
//		}
////		assert(false) : hashCycles+", "+cycleMask+", "+bitsPerCycle+", "+code;
//		return code;
//	}
	
//	/**
//	 * Generates a hash code from a kmer.
//	 * @param kmer Kmer to hash
//	 * @return Hash code
//	 */
//	private static final long hash1(long kmer){
//		long code=kmer;
//		for(int i=0, offset=0; i<hashCycles; i++, offset+=codeIncrement){
//			int x=(int)(kmer&cycleMask);
//			kmer>>=bitsPerCycle;
//			code^=codes1D[offset+x];
////			code^=Long.rotateRight(codes[i][x], (int)(code&3)); //Might be better.
//		}
////		assert(false) : hashCycles+", "+cycleMask+", "+bitsPerCycle+", "+code;
//		return code;
//	}
	
	/**
	 * Fastest version!
	 * Generates a hash code from a kmer.
	 * @param kmer Kmer to hash
	 * @return Hash code
	 */
	private static final long hash1(long kmer){
		long code=kmer;
		for(int i=0; i<codeMax; i+=codeIncrement){
			int x=(int)(kmer&cycleMask);
			kmer>>=bitsPerCycle;
			code^=codes1D[i+x];
		}
		
//		if(code>=SketchObject.minHashValue){
//			System.err.println("Made code   "+code+"\t"+Sketch.toA48(code)+"\t"+Sketch.toA48(Long.MAX_VALUE-code));
//		}else{
//			System.err.println(code+" < "+SketchObject.minHashValue);
//		}
		
		return code;
	}
	
	
//	/**
//	 * Generates a hash code from a kmer, using dual kmer lengths.
//	 * @param kmer Kmer to hash
//	 * @return Hash code
//	 */
//	private static final long hash2(final long kmer0){
//		
//		long kmer=kmer0;
//		long code=0;
//		
//		for(int i=0; i<hashCycles2; i++){
//			int x=(int)(kmer&cycleMask);
//			kmer>>=bitsPerCycle;
//			code^=codes[i][x];
//		}
//		final long code2=code^(k2mask&kmer0);
////		final long code2=(code^(k2mask&kmer0))&codeAndMask; //Distinguishes between k and k2 using the last bit.
////		if((code2&1L)==0L){return code2;} //Early exit is not as fast due to added conditional.
//		for(int i=hashCycles2; i<hashCycles; i++){
//			int x=(int)(kmer&cycleMask);
//			kmer>>=bitsPerCycle;
//			code^=codes[i][x];
//		}
////		assert(false) : hashCycles+", "+cycleMask+", "+bitsPerCycle+", "+code;
////		return code^kmer0; //For early exit version
//		return ((code&1)==1) ? code2 : code^kmer0; //Random; faster than max.
////		return ((code&1L)==0L) ? code2 : (code^kmer0)|codeOrMask; //Distinguishes between k and k2 using the last bit.
//		
//		//Functions below, based on code2, lead to underestimates of genome size.
////		return ((code2&0)==0) ? code2 : code^kmer0; //Fastest.
////		return ((code2&0)==0) ? code2 : ((code^kmer0)|codeOrMask); //Distinguishes between k and k2 using the last bit.
//	}
	
	
	/**
	 * Generates a hash code from a kmer, using dual kmer lengths.
	 * @param kmer Kmer to hash
	 * @return Hash code
	 */
	private static final long hash2(final long kmer0){
		
		long kmer=kmer0;
		long code=0;
		
		for(int i=0; i<codeMax2; i+=codeIncrement){
			int x=(int)(kmer&cycleMask);
			kmer>>=bitsPerCycle;
			code^=codes1D[i+x];
		}
		final long code2=code^(k2mask&kmer0);
		for(int i=codeMax2; i<codeMax; i+=codeIncrement){
			int x=(int)(kmer&cycleMask);
			kmer>>=bitsPerCycle;
			code^=codes1D[i+x];
		}
		return ((code&1)==1) ? code2 : code^kmer0; //Random; faster than max.
	}
	
	
//	/**
//	 * Generates a hash code from a kmer, using dual kmer lengths, allowing K2 to be a non-multiple-of-4.
//	 * Uses one additional and, xor, and lookup.
//	 * @param kmer Kmer to hash
//	 * @return Hash code
//	 */
//	private static final long hash3(final long kmer0){
//		
//		long kmer=kmer0;
//		long code=0;
//		assert(k2submask!=0);
//		
//		int i=0;
//		for(; i<hashCycles2; i++){
//			int x=(int)(kmer&cycleMask);
//			kmer>>=bitsPerCycle;
//			code^=codes[i][x];
//		}
//		final int residual=(int)(kmer&k2submask);
//		final long code2=code^(k2mask&kmer0)^codes[i][residual];
////		final long code2=(code^(k2mask&kmer0))&codeAndMask; //Distinguishes between k and k2 using the last bit.
////		if((code2&1L)==0L){return code2;} //Early exit is not as fast due to added conditional.
//		for(; i<hashCycles; i++){
//			int x=(int)(kmer&cycleMask);
//			kmer>>=bitsPerCycle;
//			code^=codes[i][x];
//		}
////		assert(false) : hashCycles+", "+cycleMask+", "+bitsPerCycle+", "+code;
////		return code^kmer0; //For early exit version
//		return ((code&1)==1) ? code2 : code^kmer0; //Random; faster than max.
////		return ((code&1L)==0L) ? code2 : (code^kmer0)|codeOrMask; //Distinguishes between k and k2 using the last bit.
//		
//		//Functions below, based on code2, lead to underestimates of genome size.
////		return ((code2&0)==0) ? code2 : code^kmer0; //Fastest.
////		return ((code2&0)==0) ? code2 : ((code^kmer0)|codeOrMask); //Distinguishes between k and k2 using the last bit.
//	}
	
	
	/**
	 * Generates a hash code from a kmer, using dual kmer lengths, allowing K2 to be a non-multiple-of-4.
	 * Uses one additional and, xor, and lookup.
	 * @param kmer Kmer to hash
	 * @return Hash code
	 */
	private static final long hash3(final long kmer0){
		
		long kmer=kmer0;
		long code=0;
		assert(k2submask!=0);
		
		int i=0;
		for(; i<codeMax2; i+=codeIncrement){
			int x=(int)(kmer&cycleMask);
			kmer>>=bitsPerCycle;
			code^=codes1D[i+x];
		}
		final int residual=(int)(kmer&k2submask);
		final long code2=code^(k2mask&kmer0)^codes1D[i+residual];
		for(; i<codeMax; i+=codeIncrement){
			int x=(int)(kmer&cycleMask);
			kmer>>=bitsPerCycle;
			code^=codes1D[i+x];
		}
		return ((code&1)==1) ? code2 : code^kmer0; //Random; faster than max.
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------      Convenience Methods     ----------------*/
	/*--------------------------------------------------------------*/
	
	public static final float aniToWkid(final double ani){
		final float wkid;
		if(ani<=0){
			wkid=0;
		}else if(k2<1){
			wkid=(float)Math.pow(ani, k);
		}else{
			wkid=0.5f*(float)(Math.pow(ani, k)+Math.pow(ani, k2));
		}
		return wkid;
	}
	
	public static final float wkidToAniExact(final double wkid, final int A, final int B, final double maxDeviation){
		assert(A>B);
		assert(maxDeviation>0);
		final double logWkid=Math.log(wkid);
		final double invA=1.0/A;
		final double invB=1.0/B;
		
		double aniUpper=Math.exp(logWkid*invA);
		double aniLower=Math.exp(logWkid*invB);
		assert(aniLower<=aniUpper);
		double aniEst=(aniLower+aniUpper)*0.5f;
		
//		if(B>16){//Fast mode
//			aniUpper=(aniUpper+3*aniEst)*0.25;
//			aniLower=(aniLower+3*aniEst)*0.25;
//		}
		
		double wkidEst=aniToWkid(aniEst, A, B);
		int iters=1;
//		System.err.println(aniLower+"\t"+aniUpper+"\t"+aniEst+"\t"+aniToWkid(aniEst, A, B)+"\t"+(wkidEst-wkid));
		while(Math.abs(wkidEst-wkid)>maxDeviation){
			if(wkidEst<wkid){aniLower=aniEst;}
			else{aniUpper=aniEst;}
			aniEst=(aniLower+aniUpper)*0.5f;
			wkidEst=aniToWkid(aniEst, A, B);
//			System.err.println(aniLower+"\t"+aniUpper+"\t"+aniEst+"\t"+wkidEst+"\t"+(wkidEst-wkid));
			iters++;
		}
//		System.err.println("Iterations: "+iters);
		return (float)aniEst;
	}
	
	public static double aniToWkid(double ani, int A, int B){
		return 0.5*(Math.pow(ani, A)+Math.pow(ani, B));
	}
	
	public static double aniToWkid(double ani, int A){
		return Math.pow(ani, A);
	}
	
	public static final float wkidToAniExact(final double wkid, final int k){
		return (float)Math.exp(Math.log(wkid)/(k));
	}
	
	public static final float wkidToAni(final double wkid){
		final float ani;
		if(wkid<=0){
			ani=0;
		}else if(k2<1){
			ani=(float)Math.exp(Math.log(wkid)/k);
		}else{
			if(EXACT_ANI){return wkidToAniExact(wkid, k, k2, 0.00005);}
			
			//This is linear interpolation which is not correct.  But for some reason it works really well.
			final double log=Math.log(wkid);
			
//			double ani1=Math.exp(log/k);
//			double ani2=Math.exp(log/k2);
//			ani=(float)(0.5*(ani1+ani2));
			
			//alternatively...  this one seems to work better.
//			ani=(float)Math.exp(2*log/(k*1.1+k2*0.9));
			ani=(float)Math.exp(2*log/(k+k2));
		}
		return ani;
	}

	public static void kmerArrayToSketchArray(long[] kmers){
		for(int i=0; i<kmers.length; i++){
			long kmer=kmers[i];
			long hash=SketchTool.hash(kmer);
			assert(hash>=minHashValue);
			hash=Long.MAX_VALUE-hash;
			kmers[i]=hash;
		}
		Shared.sort(kmers);
	}
	
	public static void hashArrayToSketchArray(long[] codes){
		for(int i=0; i<codes.length; i++){
			long hash=codes[i];
			assert(hash>=minHashValue);
			hash=Long.MAX_VALUE-hash;
			codes[i]=hash;
		}
		Shared.sort(codes);
	}
	
	public static final long genomeSizeEstimate(long min, int length) {
		if(length==0){return 0;} 
		double est=((double)Long.MAX_VALUE)*2*length/(Tools.max(min, 1));
//		if(k2>0){est*=0.5;} //This is necessary if the hash function uses max, but not random.
//		System.err.println("max="+Long.MAX_VALUE+", min="+min+", len="+length+", est="+(long)est);
//		new Exception().printStackTrace(System.err);
		return (long)Math.ceil(est);
	}
	
	public static final int toSketchSize(long genomeSizeBases, long genomeSizeKmers, long genomeSizeEstimate, int maxSketchSize){
		if(genomeSizeEstimate>0 && useSizeEstimate){
			return toSketchSizeKmers(genomeSizeEstimate, maxSketchSize);
		}
		if(genomeSizeKmers>0){
			return toSketchSizeKmers(genomeSizeKmers, maxSketchSize);
		}
		return toSketchSizeBases(genomeSizeBases, maxSketchSize);
	}
	
	private static final int toSketchSizeBases(long genomeSizeBases, int maxSketchSize) {
		return toSketchSizeKmers(Tools.max(0, genomeSizeBases-k+1), maxSketchSize);
	}
	
	private static final int toSketchSizeKmers(long genomeSizeKmers, int maxSketchSize) {
//		System.err.println(genomeSizeKmers+", "+maxSketchSize);
//		new Exception().printStackTrace();
		if(AUTOSIZE){
			if(amino){
				float linear1=Tools.min(60+0.5f*(float)Math.sqrt(genomeSizeKmers), 
						maxGenomeFractionSmall*genomeSizeKmers);
				float linear2=genomeSizeKmers*maxGenomeFraction;
				float poly=0+1f*(float)Math.sqrt(genomeSizeKmers)+90f*(float)Math.pow(genomeSizeKmers, 0.3);
				float log=Tools.max(1000, -4000+3500*(float)Math.log(Tools.max(1, genomeSizeKmers)));
				float min=Tools.min(Tools.max(linear1, linear2), poly, log);
				assert(min>=Integer.MIN_VALUE && min<=Integer.MAX_VALUE) : min;

				//			System.err.println(genomeSizeKmers+" -> "+linear1+", "+linear2+", "+poly+", "+log);

				return (int)Tools.min(genomeSizeKmers*keyFraction2, min*AUTOSIZE_FACTOR);
				//			return (int)Tools.max(1, min*AUTOSIZE_FACTOR);
			}else{
				float linear1=Tools.min(smallSketchSize+0.5f*(float)Math.sqrt(genomeSizeKmers), 
						maxGenomeFractionSmall*genomeSizeKmers-10);
				float linear2=genomeSizeKmers*maxGenomeFraction;
				float poly=0+1f*(float)Math.sqrt(genomeSizeKmers)+90f*(float)Math.pow(genomeSizeKmers, 0.3);
				float log=Tools.max(1000, -4000+3500*(float)Math.log(Tools.max(1, genomeSizeKmers))+8f*(float)Math.pow(genomeSizeKmers, 0.3));
				float min=Tools.min(Tools.max(linear1, linear2), poly, log);
				assert(min>=Integer.MIN_VALUE && min<=Integer.MAX_VALUE) : min;

				//			System.err.println(genomeSizeKmers+" -> "+linear1+", "+linear2+", "+poly+", "+log);
				
				int result=(int)Tools.min(genomeSizeKmers*keyFraction2, min*AUTOSIZE_FACTOR);
//				System.err.println(result);
				return result;
				//			return (int)Tools.max(1, min*AUTOSIZE_FACTOR);
			}
		}else{
			return Tools.min((int)(2+maxGenomeFraction*genomeSizeKmers), maxSketchSize);
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/
	
	public static final int RAW=0, HEX=1, A48=2;
	public static final char[] codingArray={'R', 'H', 'A'};
	
	static final byte[] hexTable=new byte[128];
	static {
		Arrays.fill(hexTable, (byte)-1);
		for(int i='0'; i<='9'; i++){
			hexTable[i]=(byte)(i-'0');
		}
		for(int i='A'; i<='F'; i++){
			hexTable[i]=hexTable[i+'a'-'A']=(byte)(i-'A'+10);
		}
		hexTable['x']=hexTable['X']=hexTable['-']=hexTable['+']=0;
	}
	
	public static final int ONE_SKETCH=1, PER_SEQUENCE=2, PER_TAXA=3, PER_IMG=4, PER_HEADER=5;
	
	/*--------------------------------------------------------------*/
	/*----------------      Default Locations       ----------------*/
	/*--------------------------------------------------------------*/

//	public static ArrayList<String> IMG_PATH(){return makePaths(IMG_PATH, 31);}
//	public static ArrayList<String> NT_PATH(){return makePaths(NT_PATH, 31);}
//	public static ArrayList<String> NR_PATH(){throw new RuntimeException("NR is not currently available.");}
//	public static ArrayList<String> REFSEQ_PATH(){return makePaths(REFSEQ_PATH, 31);}
//	public static ArrayList<String> SILVA_PATH(){return makePaths(SILVA_PATH, 31);}
	
	private static ArrayList<String> makePaths(String pattern, int files){
		ArrayList<String> list=new ArrayList<String>(files);
		for(int i=0; i<files; i++){
			list.add(pattern.replace("#", ""+i));
		}
		return list;
	}
	
	static final String IMG_PATH="/global/projectb/sandbox/gaag/bbtools/img/current/img#.sketch";
	static final String NT_PATH="/global/projectb/sandbox/gaag/bbtools/nt/current/taxa#.sketch";
	static final String NR_PATH="/global/projectb/sandbox/gaag/bbtools/nr/current/taxa#.sketch";
	static final String REFSEQ_PATH="/global/projectb/sandbox/gaag/bbtools/refseq/current/taxa#.sketch";
	static final String SILVA_PATH="/global/projectb/sandbox/gaag/bbtools/silva/bl_ssu_seq#.sketch";
	
	/*--------------------------------------------------------------*/
	/*----------------           Getters            ----------------*/
	/*--------------------------------------------------------------*/

	public static boolean useWhitelist() {return useWhitelist;}
	public static String blacklist() {return blacklist;}
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	public static int CODING=A48;
	public static boolean deltaOut=true;
	public static boolean rcomp=true;
	
	public static final int defaultK=31;
	public static int k=defaultK;
	public static int k2=0;
	public static int entropyK=3;
	public static boolean setK=false;
	public static boolean amino=false;
	public static boolean amino8=false;
	public static int HASH_VERSION=1;
	
	public static byte minQual=0;
	public static float minProb=0.0001f;//Q5.9; 0.00005 is Q5.6 
	
	private static final int bitsPerCycle=8; //Optimal speed for K=31 is 9bpc (15% faster than 8). 9, 10, and 11 are similar.
	private static final long cycleMask=~((-1L)<<bitsPerCycle);
	private static final long codeOrMask=1L;
	private static final long codeAndMask=~1L;
	private static final int maxCycles=(64+bitsPerCycle-1)/bitsPerCycle;
	private static final int codeIncrement=(int)(cycleMask+1);
	private static int codeMax;
	private static int codeMax2;
	
	private static long hashSeed=12345;
	private static final long[][] codes=makeCodes(maxCycles, codeIncrement, hashSeed, false);
	private static final long[] codes1D=makeCodes1D(maxCycles, codeIncrement, hashSeed, false);
	
	//These are set in postParse()
	private static int bitsPerBase=2;
	private static int basesPerCycle=bitsPerCycle/bitsPerBase;
	private static int hashCycles=64/bitsPerCycle; //Note: Needs to be variable based on k to make k2 codes compatible with k codes
	private static int hashCycles2;
	private static long k2mask;
	private static long k2submask;

//	private static Hasher hasher;
	
	public static int sketchHeapFactor=4;
//	static int minKeyOccuranceCount=1;
	
	public static int minSketchSize=3;
	public static int targetSketchSize=10000;
	public static boolean AUTOSIZE=true;
	public static float AUTOSIZE_FACTOR=1;
	public static float maxGenomeFraction=0.01f;
	public static float maxGenomeFractionSmall=0.10f;
	public static int smallSketchSize=150;
	public static boolean makeIndex=true;
	public static boolean allToAll=false;
	public static boolean compareSelf=false;
	public static int bitSetBits=2;

	private static double keyFraction=0.2;
	private static double keyFraction2=keyFraction*1.2;
	public static long minHashValue=setMinHashValue();
	public static double keyFraction(){return keyFraction;}
	public static void setKeyFraction(double d){
		assert(d>0);
		keyFraction=Tools.mid(0, d, 0.5);
		keyFraction2=Tools.mid(0, keyFraction*1.2, 0.5);
	}
	public static long setMinHashValue(){
		double mult=1-2*keyFraction;
		minHashValue=(long)(mult*Long.MAX_VALUE);
		assert(minHashValue>=0 && minHashValue<Long.MAX_VALUE);
		return minHashValue;
	}
	
	static int minFakeID=1900000000;
	static boolean hashNames=false;
	static boolean skipNonCanonical=true;
	static boolean useSizeEstimate=true;
	static boolean allowMultithreadedFastq=false;
	static boolean forceDisableMultithreadedFastq=false;
	static boolean mergePairs=false;
	
	static long sampleseed=-1L;
	
	public static TaxTree taxtree=null;
	private static String treefile=null;
	static String blacklist=null;
	static boolean useWhitelist=false;
	private static boolean postparsed=false;
	public static boolean KILL_OK=false;
	public static boolean EXACT_ANI=false; //Oddly, seems to reduce accuracy in my tests.
	
	//Needs to be last due to dependencies.
	public static DisplayParams defaultParams=new DisplayParams();
	
	public static boolean verbose2=false;
	
}
