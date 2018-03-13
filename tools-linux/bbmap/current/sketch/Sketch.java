package sketch;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import dna.AminoAcid;
import shared.Tools;
import stream.ByteBuilder;
import structures.AbstractBitSet;
import structures.IntList;
import structures.LongHashMap;
import structures.LongList;
import tax.ImgRecord;

/**
 * @author Brian Bushnell
 * @date July 7, 2016
 *
 */
public class Sketch extends SketchObject implements Comparable<Sketch> {
	
	/*--------------------------------------------------------------*/
	/*----------------         Constructors         ----------------*/
	/*--------------------------------------------------------------*/
	
	//Array should already be hashed, sorted, unique, subtracted from Long.MAX_VALUE, then reversed.
	public Sketch(long[] array_, int[] counts_){
		this(array_, counts_, -1, -1, -1, -1, -1, -1, null, null, null, null);
	}
	
	public Sketch(SketchHeap heap, boolean clearFname, boolean keepCounts){
		this(heap.toSketchArray(), null, (int)heap.taxID, heap.imgID, heap.genomeSizeBases, heap.genomeSizeKmers, heap.genomeSequences, heap.probSum, 
				heap.taxName(), heap.name0(), heap.fname(), null);
		assert(counts==null);
		if(!heap.setMode && keepCounts){
			LongHashMap map=heap.map();
			counts=new int[array.length];
			for(int i=0; i<array.length; i++){
				int count=map.get(Long.MAX_VALUE-array[i]);
				assert(count>0) : array[i]+" -> "+count+"\n"+Arrays.toString(map.values())+"\n"+Arrays.toString(map.keys());
				counts[i]=count;
			}
		}
		if(heap.setMode){heap.clearSet();}
		heap.clear(clearFname);
//		System.err.println("size="+size+", genome="+this.genomeSize+", m"); : (int)(2+maxGenomeFraction*heap.genomeSize)+", "+this.array.length;
//		assert(false) : (int)(2+maxGenomeFraction*heap.genomeSize)+", "+this.array.length;
//		assert(false) : (counts==null)+", "+heap.setMode;
	}

	public Sketch(long[] array_, int[] counts_, int taxID_, long imgID_, long gSizeBases_, long gSizeKmers_, long gSequences_, double probSum_,
			String taxName_, String name0_, String fname_, ArrayList<String> meta_){
		array=array_;
		counts=counts_;
		assert(counts==null || array==null || counts.length==array.length) : (array==null ? "null" : array.length)+", "+(counts==null ? "null" : counts.length);
		taxID=taxID_;
		imgID=imgID_;
		genomeSizeBases=gSizeBases_;
		genomeSizeKmers=gSizeKmers_;
		genomeSequences=gSequences_;
		probCorrect=probSum_<=0 ? 0f : (float)(probSum_/(Tools.max(genomeSizeKmers, 1)));
		
		taxName=fix(taxName_);
		name0=fix(name0_);
		fname=fix(fname_);
		meta=fix(meta_);
		
		if(ImgRecord.imgMap!=null && imgID>=0 && taxID<0){
			ImgRecord record=ImgRecord.imgMap.get(imgID);
			if(record!=null){
				if(record.name!=null && taxName==null){taxName=record.name;}
				taxID=record.taxID;
			}
		}
	}
	
	void addMeta(String s){
		s=fixMeta(s);
		if(s==null){return;}
		if(meta==null){meta=new ArrayList<String>(1);}
		meta.add(s);
	}
	
	void setMeta(ArrayList<String> list){
		assert(meta==null);
		meta=fix(list);
	}
	
	private static String fix(String s){
		if(s==null){return null;}
		return s.replace('\t', ' ');
	}
	
	private static String fixMeta(String s){
		if(s==null){return null;}
		int colon=s.indexOf(':');
		assert(colon>=0);
		if(s.length()==colon+5 && s.endsWith(":null")){return null;}
		return fix(s);
	}
	
	private static ArrayList<String> fix(ArrayList<String> list){
		if(list==null || list.isEmpty()){return null;}
		for(int i=0; i<list.size(); i++){
			String s=list.get(i);
			s=fixMeta(s);
			if(s==null){
				list.remove(i);
				i--;
			}else{
				list.set(i, s);
			}
		}
		if(list==null || list.isEmpty()){return null;}
		list.trimToSize();
		Collections.sort(list);
		return list;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	public void add(Sketch other, int maxlen){
		final long[] a=array;
		final long[] b=other.array;
		if(maxlen<1){
			assert(false);
			maxlen=1000000;
		}
		LongList list=new LongList(Tools.min(maxlen, a.length+b.length));
		
		for(int i=0, j=0; i<a.length && j<b.length; ){
			final long ka=a[i], kb=b[j];
			if(ka==kb){//match
				list.add(ka);
				i++;
				j++;
			}else if(ka<kb){
				list.add(ka);
				i++;
			}else{
				list.add(kb);
				j++;
			}
			if(list.size()>=maxlen){break;}
		}
		
		if(array.length==list.size()){
			for(int i=0; i<list.size; i++){
				array[i]=list.array[i];
			}
		}else{
			array=list.toArray();
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Set Operations        ----------------*/
	/*--------------------------------------------------------------*/
	
	public static final Sketch intersection(Sketch sa, Sketch sb){
		Sketch shared=intersection(sa.array, sb.array, sa.counts);
		if(shared!=null){
			shared.taxID=sb.taxID;
			shared.taxName=sb.taxName;
			shared.name0=sb.name0;
			shared.fname=sb.fname;
			shared.meta=sb.meta;
			shared.imgID=sb.imgID;
			shared.spid=sb.spid;
		}
		return shared;
	}
	
	public static final Sketch intersection(long[] a, long[] b, int[] aCounts){
		int i=0, j=0, matches=0;
		LongList ll=new LongList();
		IntList il=new IntList();
		for(; i<a.length && j<b.length; ){
			final long ka=a[i], kb=b[j];
			if(ka==kb){
				matches++;
				ll.add(ka);
				if(aCounts!=null){
					il.add(aCounts[i]);
				}
				i++;
				j++;
			}else if(ka<kb){
				i++;
			}else{
				j++;
			}
		}
		if(matches<1){return null;}
			
		return new Sketch(ll.toArray(), il.size>0 ? il.toArray() : null);
	}
	
	public static final Sketch union(Sketch sa, Sketch sb){
		Sketch shared=union(sa.array, sb.array, sa.counts, sb.counts);
		if(shared!=null){
			shared.taxID=sa.taxID;
			shared.taxName=sa.taxName;
			shared.name0=sa.name0;
			shared.fname=sa.fname;
			shared.meta=sa.meta;
			shared.imgID=sa.imgID;
			shared.spid=sa.spid;
		}
		return shared;
	}
	
	public static final Sketch union(long[] a, long[] b, int[] aCounts, int[] bCounts){
		int i=0, j=0, matches=0;
		LongList ll=new LongList();
		IntList il=new IntList();
		for(; i<a.length && j<b.length; ){
			final long ka=a[i], kb=b[j];
			if(ka==kb){
				matches++;
				ll.add(ka);
				if(aCounts!=null && bCounts!=null){
					il.add(aCounts[i]+bCounts[i]);
				}
				i++;
				j++;
			}else if(ka<kb){
				ll.add(ka);
				if(aCounts!=null && bCounts!=null){
					il.add(aCounts[i]);
				}
				i++;
			}else{
				ll.add(kb);
				if(aCounts!=null && bCounts!=null){
					il.add(bCounts[i]);
				}
				j++;
			}
		}
		if(matches<1){return null;}
			
		return new Sketch(ll.toArray(), il.size>0 ? il.toArray() : null);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Comparison          ----------------*/
	/*--------------------------------------------------------------*/
	
	public int countMatches(Sketch other, CompareBuffer buffer, AbstractBitSet present, boolean fillPresent, int[][] taxHits, int contamLevel){
		return countMatches(array, other.array, counts, other.counts, other.taxID, buffer, present, fillPresent, taxHits, contamLevel);
	}
	
	public static final int countMatches(long[] a, long[] b, int[] aCounts, int[] bCounts, int bid, 
			CompareBuffer buffer, AbstractBitSet present, boolean fillPresent, int[][] taxHits, int contamLevel){

//		assert(fillPresent) : bid+", "+minFakeID+", "+(taxHits!=null);
		
		if(bid>0 && bid<minFakeID && taxHits!=null){
			bid=taxtree.getIdAtLevelExtended(bid, contamLevel);
		}else{
			bid=-1;
		}
		
//		assert(false) : (buffer==null)+", "+fillPresent+", "+present.cardinality();
		assert(a.length>0 && b.length>0);
		
		//Kmers hitting this reference
		int matches=0;
		
		//Kmers hitting this reference and others
		int multiMatches=0;
		
		//Kmers hitting nothing
		int noHits=0;
		
		//Kmers hitting some organism but not this reference
		int contamHits=0;
		
		//Kmers hitting something in this taxa, but not this reference 
		int sameTax=0;
		
		//Kmers hitting this organism and no other taxa
		int unique2=0;
		
		//Kmers hitting only this taxa but not this organism (this count may not include everything due to early exit) 
		int unique3_temp=0;
		
		//Kmers hitting multiple organisms but not this reference
		int multiContamHits=0;
		
		//Sum of query counts for shared kmers
		long depthSum=0;
		
		//Sum of query counts for shared kmers divided by ref counts for those kmers
		double depthSum2=0;//Slow, but necessary.
		
		int i=0, j=0;
		assert(present==null || present.capacity()==a.length);
//		assert(false) : buffer.rbs.capacity()+", "+buffer.rbs+", "+present;
		if(present!=null){
			if(fillPresent){
				for(; i<a.length && j<b.length; ){
					final long ka=a[i], kb=b[j];
					if(ka==kb){
						present.increment(i);
						matches++;
						if(aCounts!=null){
							depthSum+=aCounts[i];
							if(bCounts!=null){
								depthSum2+=aCounts[i]/(double)bCounts[j];
							}
						}
						i++;
						j++;
					}else if(ka<kb){
						i++;
					}else{
						j++;
					}
				}
			}else{
				for(; i<a.length && j<b.length; ){
					final long ka=a[i], kb=b[j];
					if(ka==kb){
						final int count=present.getCount(i);
						if(count>1){
							multiMatches++;
						}
						
						matches++;
						if(aCounts!=null){
							depthSum+=aCounts[i];
							if(bCounts!=null){
								depthSum2+=aCounts[i]/(double)bCounts[j];
							}
						}
						if(bid>0){
							int[] taxHitsRow=taxHits[i];
							if(taxHitsRow!=null && taxHitsRow.length==1 && taxHitsRow[0]==bid){unique2++;}
						}
						
						i++;
						j++;
					}else if(ka<kb){
						final int count=present.getCount(i);
						if(count>0){
							contamHits++;
							if(count>1){
								multiContamHits++;
							}
						}else{
							noHits++;
						}
						
						if(bid>0){
							int[] taxHitsRow=taxHits[i];
							if(taxHitsRow!=null){
								if(taxHitsRow!=null && taxHitsRow.length==1 && taxHitsRow[0]==bid){unique3_temp++;}
								for(int tid : taxHitsRow){
									if(tid==bid){
										sameTax++;
										break;
									}
								}
							}
						}
						
						i++;
					}else{
						j++;
					}
				}
				
				//For the remaining query kmers, we don't know whether the reference sketch would have shared them had it been longer.
				//This section can be disabled to prevent them from being displayed.
				if(bid>0 && i<a.length-1){
					for(; i<a.length; i++){
						int[] taxHitsRow=taxHits[i];
						if(taxHitsRow!=null){
							if(taxHitsRow!=null && taxHitsRow.length==1 && taxHitsRow[0]==bid){unique3_temp++;}
						}
					}
				}
			}
		}else{
			for(; i<a.length && j<b.length; ){
				final long ka=a[i], kb=b[j];
				if(ka==kb){
					matches++;
					if(aCounts!=null){depthSum+=aCounts[i];}
					i++;
					j++;
				}else if(ka<kb){
					i++;
				}else{
					j++;
				}
			}
		}
		
//		if(taxHits!=null){
//			System.err.println("matches="+matches+", noHits="+noHits+", contamHits="+contamHits+", sameTax="+sameTax+", multiContamHits="+multiContamHits);
//		}

//		assert(bid<1 || unique2>=(matches-multiMatches)) : bid+", "+unique2+", "+unique3_temp+", "+matches+", "+multiMatches;
//		assert(matches<1000 || multiMatches==0) : bid+", "+unique2+", "+unique3_temp+", "+matches+", "+multiMatches+", "+fillPresent;
		
		if(buffer!=null){
			buffer.set(matches, multiMatches, unique2, unique2+unique3_temp, noHits, contamHits, contamHits-sameTax, multiContamHits, i, j, a.length, b.length, depthSum, depthSum2);
		}
		return matches;
	}
	
//	public float identity(Sketch b, float[] ret){
//		if(ret!=null){Arrays.fill(ret, 0);}
//		return identityWeighted(array, b.array, ret);
//	}
//	
//	public static float identity(long[] a, long[] b){
//		int matches=countMatches(a, b);
//		return matches/(float)(Tools.max(1, Tools.min(a.length, b.length)));
//	}
	
	@Override
	public int hashCode(){
		long gSize=genomeSizeKmers>0 ? genomeSizeKmers : genomeSizeBases;
		int code=(int) ((gSize^taxID^imgID^(name0==null ? 0 : name0.hashCode()))&Integer.MAX_VALUE);
//		System.err.println(code+", "+gSize+", "+taxID+", "+imgID+", "+name0);
		return code;
	}
	
	@Override
	public int compareTo(Sketch b){
		if(this==b){return 0;}
		if(taxID>-1 && b.taxID>-1){return taxID-b.taxID;}
		int x=taxName.compareTo(b.taxName);
		if(x!=0){return x;}
		if(name0!=null && b.name0!=null){return name0.compareTo(b.name0);}
		return name0!=null ? 1 : b.name0!=null ? -1 : 0;
	}
	
	@Override
	public boolean equals(Object b){
		if(this==b){return true;}
		if(b==null || this.getClass()!=b.getClass()){return false;}
		return equals((Sketch)b);
	}
	
	public boolean equals(Sketch b){
		return compareTo(b)==0;
	}
	
	public ByteBuilder toHeader(){
		ByteBuilder sb=new ByteBuilder();
		return toHeader(sb);
	}
	
	public ByteBuilder toHeader(ByteBuilder sb){
		sb.append("#SZ:").append(array.length);
		sb.append("\tCD:");
		sb.append(codingArray[CODING]);
		if(deltaOut){sb.append('D');}
		if(counts!=null){sb.append('C');}
		if(amino){sb.append('M');}
		if(amino8){sb.append('8');}

		if(k!=defaultK || k2!=0){
			sb.append("\tK:").append(k);
			if(k2!=0){sb.append(",").append(k2);}
		}
		if(HASH_VERSION>1){sb.append("H:").append(HASH_VERSION);}
		
		if(genomeSizeBases>0){sb.append("\tGS:").append(genomeSizeBases);}
		if(genomeSizeKmers>0){sb.append("\tGK:").append(genomeSizeKmers);}
		final long ge=genomeSizeEstimate();
		if(ge>0){sb.append("\tGE:").append(ge);}
		if(genomeSequences>0){sb.append("\tGQ:"+genomeSequences);}
		if(probCorrect>0){sb.append("\tPC:"+String.format("%.4f", probCorrect));}
		if(taxID>=0){sb.append("\tID:").append(taxID);}
		if(imgID>=0){sb.append("\tIMG:").append(imgID);}
		if(spid>0){sb.append("\tSPID:").append(spid);}
		if(fname!=null){sb.append("\tFN:").append(fname);}
		if(taxName!=null){sb.append("\tNM:").append(taxName);}
		if(name0!=null){sb.append("\tNM0:").append(name0);}
		if(meta!=null){
			for(String s : meta){
				sb.append("\tMT_").append(s);
			}
		}
		return sb;
	}
	
	public ByteBuilder toBytes(){
		return toBytes(new ByteBuilder());
	}
	
	public ByteBuilder toBytes(ByteBuilder sb){
		if(CODING==A48 && deltaOut){return toBytesA48D(sb);}
		long prev=0;
		toHeader(sb);
		sb.append("\n");
		byte[] temp=null;
		if(CODING==A48){temp=new byte[12];}
		for(int i=0; i<array.length; i++){
			long key=array[i];
			int count=(counts==null ? 1 : counts[i]);
			long x=key-prev;
			if(CODING==A48){
				appendA48(x, sb, temp);
				if(count>1){
					sb.append('\t');
					appendA48(count-1, sb, temp);
				}
				sb.append('\n');
			}else if(CODING==HEX){
				sb.append(Long.toHexString(x)).append('\n');
			}else if(CODING==RAW){
				sb.append(x).append('\n');
			}else{
				assert(false);
			}
			if(deltaOut){prev=key;}
		}
		return sb;
	}
	
	//This is to make the common case fast
	private ByteBuilder toBytesA48D(ByteBuilder sb){
		assert(CODING==A48 && deltaOut);
		long prev=0;
		toHeader(sb);
		sb.append("\n");
		final byte[] temp=new byte[12];

		if(counts==null){
			for(int i=0; i<array.length; i++){
				long key=array[i];
				long x=key-prev;
				if(CODING==A48){
					appendA48(x, sb, temp);
					sb.append('\n');
				}
				prev=key;
			}
		}else{
			for(int i=0; i<array.length; i++){
				long key=array[i];
				int count=counts[i];
				long x=key-prev;
				if(CODING==A48){
					appendA48(x, sb, temp);
					if(count>1){
						sb.append('\t');
						appendA48(count-1, sb, temp);
					}
					sb.append('\n');
				}
				prev=key;
			}
		}
		return sb;
	}
	
	public static final void appendA48(long value, ByteBuilder sb, byte[] temp){
		int i=0;
//		long value=value0;
		while(value!=0){
			byte b=(byte)(value&0x3F);
//			assert(i<temp.length) : i+", "+temp.length+", "+value0;
			temp[i]=b;
			value=value>>6;
			i++;
		}
		if(i==0){
			sb.append((byte)'0');
		}else{
			for(i--;i>=0;i--){
				sb.append((char)(temp[i]+48));
			}
		}
	}
	
	public static final String toA48(long value){
		int i=0;
//		long value=value0;
		StringBuilder sb=new StringBuilder(12);
		while(value!=0){
			byte b=(byte)(value&0x3F);
//			assert(i<temp.length) : i+", "+temp.length+", "+value0;
			sb.append((char)(b+48));
			value=value>>6;
			i++;
		}
		if(i==0){
			sb.append((byte)'0');
		}else{
			sb.reverse();
		}
		return sb.toString();
	}
	
	public String toString(){
		return toBytes().toString();
	}
	
	public static long parseA48(String line){
		if(line.length()==0){return 0;}
		long x=0;
		for(int i=0; i<line.length(); i++){
			x<<=6;
			long c=line.charAt(i);
			x|=(c-48);
		}
		return x;
	}
	
	/** Parses coverage too */
	public static long parseA48C(String line, IntList covList){
		if(line.length()==0){
			covList.add(1);
			return 0;
		}
		long key=0, cov=0;
		int i=0, len=line.length();
		for(; i<len; i++){
			long c=line.charAt(i);
			if(c<48){break;}
			key<<=6;
			key|=(c-48);
		}
		for(i++; i<len; i++){
			long c=line.charAt(i);
			cov<<=6;
			cov|=(c-48);
		}
		covList.add((int)(cov+1));
		return key;
	}
	
	public static long parseHex(String line){
		if(line.length()==0){return 0;}
		long x=0;
		for(int i=0; i<line.length(); i++){
			x<<=4;
			x|=hexTable[line.charAt(i)];
		}
		if(line.charAt(0)=='-'){x*=-1;}
		return x;
	}
	
	public static long parseA48(byte[] line){
		if(line.length==0){return 0;}
		long x=0;
		for(byte b : line){
			x<<=6;
			x|=(((long)b)-48);
		}
		return x;
	}
	
	public static long parseNuc(String line){
		return parseNuc(line.getBytes());
	}
	
	/** Returns the maximal key in the sequence */
	public static long parseNuc(byte[] bases){
		if(bases.length<k){return -1;}
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift)); //Conditional allows K=32
		
		long kmer=0, rkmer=0;
		int len=0;
		
		long key=-1;
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			long x=AminoAcid.baseToNumber[b];
			long x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=(rkmer>>>2)|(x2<<shift2);
			if(x<0){len=0;}else{len++;}
			if(len>=k){
				long z=Tools.max(kmer, rkmer);
				long hash=SketchTool.hash(z);
				key=Tools.max(key, hash);
			}
		}
		return key<minHashValue ? -1 : Long.MAX_VALUE-key;
	}
	
	/** Parses coverage too */
	public static long parseA48C(byte[] line, IntList covList){
		if(line.length==0){
			covList.add(1);
			return 0;
		}
		long key=0, cov=0;
		int i=0, len=line.length;
		for(; i<len; i++){
			long b=line[i];
			if(b<48){break;}
			key<<=6;
			key|=(b-48);
		}
		for(i++; i<len; i++){
			long b=line[i];
			cov<<=6;
			cov|=(b-48);
		}
		covList.add((int)(cov+1));
		return key;
	}
	
	public static long parseHex(byte[] line){
		if(line.length==0){return 0;}
		long x=0;
		for(byte b : line){
			x<<=4;
			x|=hexTable[b];
		}
		if(line[0]=='-'){x*=-1;}
		return x;
	}
	
	public long genomeSizeEstimate() {
		return array.length==0 ? 0 : Tools.min(genomeSizeKmers, genomeSizeEstimate(array[array.length-1], array.length));
	}
	
	public String name(){return taxName!=null ? taxName : name0!=null ? name0 : fname;}
	public String taxName(){return taxName;}
	public String name0(){return name0;}
	public String fname(){return fname;}
	public int length(){return array.length;}
	public void setTaxName(String s){taxName=s;}
	public void setName0(String s){name0=s;}
	public void setFname(String s){
//		assert(!s.endsWith("sketch")) : s; //123
		fname=s;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	public void makeBitSets(boolean printContam, boolean index){
		assert(compareBitSet==null && indexBitSet==null);
		if(!printContam){return;}
		compareBitSet=AbstractBitSet.make(length(), bitSetBits);
		if(index){indexBitSet=AbstractBitSet.make(length(), bitSetBits);}
	}
	
	public void addToBitSet(AbstractBitSet rbs){
		compareBitSet.add(rbs);
	}
	
	public AbstractBitSet compareBitSet(){return compareBitSet;}
	
	
	public AbstractBitSet indexBitSet(){return indexBitSet;}
	
	public void mergeBitSets(){
		assert(!mergedBitSets);
		if(compareBitSet!=null && indexBitSet!=null){
			compareBitSet.setToMax(indexBitSet);
		}
		indexBitSet=null;
		mergedBitSets=true;
	}
	
	public boolean merged(){return mergedBitSets;}
	
	public long[] array;
	int[] counts;
	public int taxID;
	public final long genomeSequences;
	public final long genomeSizeBases;
	public final long genomeSizeKmers;
	public final float probCorrect;
	private String taxName;
	private String name0;
	private String fname;
	ArrayList<String> meta;
	
	//TODO: These should move to SketchResults.
	private AbstractBitSet compareBitSet; //Used for comparison
	private AbstractBitSet indexBitSet;
	
	//Extended information
	public long imgID=-1;
	public long spid=-1;
//	public String seqUnitName=null;
	
	private boolean mergedBitSets=false; //Temporary for debugging
	
	
	
}
