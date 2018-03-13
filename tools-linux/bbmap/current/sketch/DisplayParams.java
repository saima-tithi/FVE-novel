package sketch;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.Locale;

import shared.Colors;
import shared.Tools;
import tax.PrintTaxonomy;
import tax.TaxFilter;
import tax.TaxNode;
import tax.TaxTree;

public class DisplayParams implements Cloneable {
	
	@Override
	public DisplayParams clone(){
		try {
			DisplayParams copy=(DisplayParams)super.clone();
			if(taxFilter!=null){
				copy.taxFilter=taxFilter.deepCopy();
			}
			return copy;
		} catch (CloneNotSupportedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			throw new RuntimeException();
		}
	}
	
	public DisplayParams parseDoubleHeader(String s){
		if(!s.startsWith("##")){return this;}
		StringBuilder sb=new StringBuilder();
		for(int i=2; i<s.length(); i++){
			char c=s.charAt(i);
			if(c=='\n'){break;}
			sb.append(c);
		}
		return parseDoubleHeaderLine(sb.toString());
	}
	
	public DisplayParams parseDoubleHeaderLine(String line) {
		if(line.startsWith("##")){line=line.substring(2);}
		else{assert(!line.startsWith("#")) : line;}
		if(line.length()<1){return this;}
		
		DisplayParams params=(DisplayParams) this.clone();
		
		String[] args=line.split(" ");
		for(String arg : args){
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b==null || b.equalsIgnoreCase("null")){b=null;}
			while(a.startsWith("-")){a=a.substring(1);} //Strip leading hyphens
			
			boolean x=params.parse(arg, a, b);
//			assert(x) : "Unknown parameter "+arg+"\n"+line;
			if(!x){System.err.println("Warning: Unknown parameter "+arg);}
		}
		if(SketchObject.verbose2){System.err.println("Made it to post-parse.  taxFilter="+params.taxFilter);}
		params.postParse(true);
		if(SketchObject.verbose2){System.err.println("Passed post-parse.  taxFilter="+params.taxFilter);}
		
		return params;
	}
	
	public boolean parse(String arg, String a, String b){
	
		if(a.equals("minhits")  || a.equals("hits")){
			minHits=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("minwkid") || a.equalsIgnoreCase("wkid")){
			minWKID=Float.parseFloat(b);
			if(minWKID>1){minWKID/=100;}
			assert(minWKID<=1) : "minWKID should between 0 and 1";
		}else if(a.equalsIgnoreCase("minid") || a.equalsIgnoreCase("id") || a.equalsIgnoreCase("minani") || a.equalsIgnoreCase("ani")){
			minANI=Float.parseFloat(b);
			if(minANI>1){minANI/=100;}
			assert(minANI<=1) : "minANI should between 0 and 1";
			if(minANI>0){
				minWKID=(float)Tools.max(minWKID, Comparison.aniToWkid(minANI, 32));//Lowest possible minWKID for this ANI
			}
		}else if(a.equals("records") || a.equals("maxrecords") || a.equals("results")){
			maxRecords=Integer.parseInt(b);
			assert(maxRecords>=1) : "Max records must be at least 1.";
		}else if(a.equals("format")){
			format=Integer.parseInt(b);
		}else if(a.equals("level") || a.equals("taxlevel") || a.equals("minlevel")){
			taxLevel=TaxTree.parseLevel(b);//TODO: Change to extended
		}
		
		else if(a.equalsIgnoreCase("printtax") || a.equalsIgnoreCase("printtaxa")){
			printTax=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printoriginalname") || a.equalsIgnoreCase("printseqname") || a.equalsIgnoreCase("printname0") || a.equals("pn0")){
			printOriginalName=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printfilename") || a.equalsIgnoreCase("printfname")){
			printFileName=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printimg")){
			printImg=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printcompleteness") || a.equalsIgnoreCase("completeness") || a.equalsIgnoreCase("printcomplt")){
			printCompleteness=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printani") || a.equalsIgnoreCase("ani")){
			printAni=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printscore") || a.equalsIgnoreCase("score")){
			printScore=Tools.parseBoolean(b);
		}
		
		else if(a.equalsIgnoreCase("trackcounts")){
			trackCounts=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printdepth") || a.equalsIgnoreCase("depth")){
			printDepth=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printdepth2") || a.equalsIgnoreCase("depth2")){
			printDepth2=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("actualdepth")){
			printActualDepth=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printvolume") || a.equalsIgnoreCase("volume")){
			printVolume=Tools.parseBoolean(b);
		}
		
		else if(a.equalsIgnoreCase("sortByDepth")){
			boolean x=Tools.parseBoolean(b);
			if(x){comparator=Comparison.depthComparator;}
		}else if(a.equalsIgnoreCase("sortByDepth2")){
			boolean x=Tools.parseBoolean(b);
			if(x){comparator=Comparison.depth2Comparator;}
		}else if(a.equalsIgnoreCase("sortByVolume")){
			boolean x=Tools.parseBoolean(b);
			if(x){comparator=Comparison.volumeComparator;}
		}else if(a.equalsIgnoreCase("sortByScore")){
			boolean x=Tools.parseBoolean(b);
			if(x){comparator=Comparison.scoreComparator;}
		}
		
		else if(a.equalsIgnoreCase("printUMatches") || a.equalsIgnoreCase("printUHits") || a.equalsIgnoreCase("printUnique")){
			printUnique=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printUMatches2") || a.equalsIgnoreCase("printUnique2") || a.equalsIgnoreCase("unique2")){
			printUnique2=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printUMatches3") || a.equalsIgnoreCase("printUnique3") || a.equalsIgnoreCase("unique3")){
			printUnique3=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printUContam")){
			printUContam=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printNoHit")){
			printNoHit=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("contamhits") || a.equalsIgnoreCase("contam") || a.equalsIgnoreCase("printcontam")){
			printContam=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("contamhits2") || a.equalsIgnoreCase("contam2") || a.equalsIgnoreCase("printcontam2")){
			if(b==null || b.length()<1){
				printContam2=true;
			}else if(Character.isDigit(b.charAt(0)) || b.charAt(0)=='-'){
				contamLevel=Tools.max(0, TaxTree.levelToExtended(Integer.parseInt(b)));
				printContam2=true;
			}else if(TaxTree.levelMapExtendedContains(b)){
				contamLevel=TaxTree.stringToLevelExtended(b);
				printContam2=true;
			}else{
				printContam2=Tools.parseBoolean(b);
			}
		}else if(a.equalsIgnoreCase("contamLevel")){
			if(Character.isDigit(b.charAt(0)) || b.charAt(0)=='-'){
				contamLevel=Tools.max(0, TaxTree.levelToExtended(Integer.parseInt(b)));
				printContam2=true;
			}else if(TaxTree.levelMapExtendedContains(b)){
				contamLevel=TaxTree.stringToLevelExtended(b);
				printContam2=true;
			}
		}
		
		else if(a.equalsIgnoreCase("printMatches")){
			printMatches=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printLength")){
			printLength=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printTaxID")){
			printTaxID=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printGSize")){
			printGSize=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printGKmers")){
			printGKmers=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printTaxName")){
			printTaxName=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printGSeqs")){
			printGSeqs=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printGBases")){
			printGBases=Tools.parseBoolean(b);
		}
		
		else if(a.equalsIgnoreCase("minEntropy") || a.equalsIgnoreCase("entropy") || a.equalsIgnoreCase("efilter")){
			minEntropy=Float.parseFloat(b);
		}
		
		else if(a.equalsIgnoreCase("printColors") || a.equalsIgnoreCase("colors") || a.equalsIgnoreCase("color")){
//			System.err.println("Parsing '"+arg+"'"); //123
			if(b==null || b.length()<1){
				printColors=true;
			}else if(b.equalsIgnoreCase("t") || b.equalsIgnoreCase("true")){
				printColors=true;
			}else if(b.equalsIgnoreCase("f") || b.equalsIgnoreCase("false")){
				printColors=false;
			}else{
				printColors=true;
				if(Character.isDigit(b.charAt(0)) || b.charAt(0)=='-'){
					colorLevel=Tools.max(0, TaxTree.levelToExtended(Integer.parseInt(b)));
				}else{
					colorLevel=TaxTree.stringToLevelExtended(b);
				}
			}
			setColors=true;
//			System.err.println("Parsed "+arg); //123
		}else if(a.equalsIgnoreCase("colorLevel")){
//			System.err.println("Parsing '"+arg+"'"); //123
			if(Character.isDigit(b.charAt(0)) || b.charAt(0)=='-'){
				colorLevel=Tools.max(0, TaxTree.levelToExtended(Integer.parseInt(b)));
			}else{
				colorLevel=TaxTree.stringToLevelExtended(b);
			}
//			System.err.println("Parsed "+arg); //123
		}
		
		else if(a.equalsIgnoreCase("printRefDivisor") || a.equalsIgnoreCase("printRDiv")){
			printRefDivisor=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printQueryDivisor") || a.equalsIgnoreCase("printQDiv")){
			printQueryDivisor=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printRefSize") || a.equalsIgnoreCase("printRSize")){
			printRefSize=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printQuerySize") || a.equalsIgnoreCase("printQSize")){
			printQuerySize=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printContamHits") || a.equalsIgnoreCase("printCHits")){
			printContamHits=Tools.parseBoolean(b);
		}
		
		else if(a.equalsIgnoreCase("printIntersection") || a.equalsIgnoreCase("intersection") || a.equalsIgnoreCase("intersect")){
			printIntersection=Tools.parseBoolean(b);
		}
		
		else if(a.equalsIgnoreCase("printAll")){
			if(Tools.parseBoolean(b)){
				setPrintAll();
			}
		}
		
		else if(a.equals("samplerate")){
			samplerate=Float.parseFloat(b);
		}else if(a.equals("reads")){
			reads=Tools.parseKMG(b);
		}else if(a.equals("mode") || a.equalsIgnoreCase("single") || a.equalsIgnoreCase("singlesketch") || a.equalsIgnoreCase("onesketch") 
				|| a.equalsIgnoreCase("persequence") || a.equalsIgnoreCase("pertaxa") || a.equalsIgnoreCase("perheader")){
			mode=SketchObject.parseMode(arg, a, b);
		}
		
		//For format 3
		else if(a.equalsIgnoreCase("useTaxidName") || a.equalsIgnoreCase("useTaxidAsName")){
			useTaxidName=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("useImgName") || a.equalsIgnoreCase("useImgAsName")){
			useImgName=Tools.parseBoolean(b);
		}else if(a.equalsIgnoreCase("useTaxName") || a.equalsIgnoreCase("useTaxAsName")){
			useTaxName=Tools.parseBoolean(b);
		}
		
		else if(a.equalsIgnoreCase("taxfilter") || a.equalsIgnoreCase("taxfilterset")){
			if(b==null){taxFilter=null;}
			else{
				if(taxFilter==null){taxFilter=new TaxFilter(SketchObject.taxtree);}
				taxFilter.clearSet();
				taxFilter.makeSet();
				taxFilter.addNumbers(b, false);
//				System.err.println("A:\t"+this);
			}
		}else if(a.equalsIgnoreCase("taxfilterlevel")){//TODO:  Change to extended
			int temp=TaxTree.parseLevel(b);
			if(taxFilter==null){taxFilter=new TaxFilter(SketchObject.taxtree);}
			taxFilter.setLevel(temp, false);
//			System.err.println("B:\t"+this);
		}else if(a.equalsIgnoreCase("taxfilterinclude") || a.equalsIgnoreCase("taxfilterexclude")){
			boolean temp=Tools.parseBoolean(b);
			if(a.equalsIgnoreCase("taxfilterexclude")){temp=!temp;}
			if(taxFilter==null){taxFilter=new TaxFilter(SketchObject.taxtree);}
			taxFilter.setInclude(temp);
//			System.err.println("C:\t"+this);
		}else if(a.equalsIgnoreCase("taxfilterstring")){
			if(taxFilter==null){taxFilter=new TaxFilter(SketchObject.taxtree);}
			taxFilter.setContainsString(b);
//			System.err.println("D:\t"+this);
		}
		
		else if(a.equalsIgnoreCase("minkmercount") || a.equalsIgnoreCase("minkeycount") || a.equalsIgnoreCase("mincount") || a.equalsIgnoreCase("minKeyOccuranceCount")){
			minKeyOccuranceCount=Integer.parseInt(b);
		}
		
		//TODO:  Eventually remove support for "amino" and "k" and just support "hamino" and "hk"
		
		//Parameters for compatibility verification
		else if(a.equalsIgnoreCase("k") || a.equalsIgnoreCase("hk")){
//			System.err.println("A: k="+k+", k2="+k2+", arg="+arg);
			if(b.indexOf(',')>=0){
				String[] split=b.split(",");
				assert(split.length==2) : "\nBad argument "+arg+"\n"+b+"\n";
				int x=Integer.parseInt(split[0]);
				int y=Integer.parseInt(split[1]);
				k=Tools.max(x, y);
				k2=Tools.min(x, y);
//				System.err.println("B: k="+k+", k2="+k2+", split="+Arrays.toString(split));
			}else{
				k=Integer.parseInt(b);
//				System.err.println("C: k="+k+", k2="+k2);
			}
		}else if(a.equals("hashversion") || a.equals("hv")){
			hashVersion=Integer.parseInt(b);
		}else if(a.equals("amino") || a.equals("hamino")){
			amino=Tools.parseBoolean(b);
		}
		
		else{
			return false;
		}
		return true;
	}
	
	public void postParse(boolean requireTree){
		assert(!postParsed);
		synchronized(this){
			if(postParsed){return;}
			if(taxFilter!=null){
				if(taxFilter.size()==0 && taxFilter.containsString()==null){
					System.err.println("Eliminating empty TaxFilter.");
					taxFilter=null;
				}
			}
			if(taxFilter!=null && requireTree){
				assert(SketchObject.taxtree!=null) : "No taxtree loaded.";
				taxFilter.setTree(SketchObject.taxtree);
				taxFilter.promote();
			}
			postParsed=true;
		}
	}
	
	public String toString(){
		StringBuilder sb=new StringBuilder();
		sb.append("##");
		sb.append("hits=").append(minHits);
		sb.append(" wkid=").append(String.format(Locale.ROOT, "%.5f",minWKID));
		if(minANI>0){sb.append(" id=").append(String.format(Locale.ROOT, "%.5f",minANI));}
		sb.append(" records=").append(maxRecords);
		sb.append(" format=").append(format);
		sb.append(" level=").append(taxLevel);
		
		if(k!=SketchObject.defaultK || k2!=0 || k!=SketchObject.k || k2!=SketchObject.k2){
			assert(k>0 && k2>=0 && k2<k) : "Bad values for k: "+k+", "+k2+", "+SketchObject.k+", "+SketchObject.k2;
			assert(SketchObject.k>0 && SketchObject.k2>=0 && SketchObject.k2<SketchObject.k) : "Bad values for k: "+k+", "+k2+", "+SketchObject.k+", "+SketchObject.k2;
			sb.append(" hk=").append(SketchObject.k).append(',').append(SketchObject.k2);
		}
		if(SketchObject.amino){sb.append(" hamino=").append(SketchObject.amino);} //TODO: This conflicts with Parser flag
		if(SketchObject.HASH_VERSION>1){sb.append(" hashversion=").append(SketchObject.HASH_VERSION);}

		if(true || printTax!=default_printTax){sb.append(" printTax=").append(printTax);}
		if(true || printOriginalName!=default_printOriginalName){sb.append(" pn0=").append(printOriginalName);}
		if(true || printFileName!=default_printFileName){sb.append(" printfname=").append(printFileName);}
		if(true || printImg!=default_printImg){sb.append(" printImg=").append(printImg);}
		if(true || printAni!=default_printAni){sb.append(" printAni=").append(printAni);}
		if(true || printCompleteness!=default_printCompleteness){sb.append(" printCompleteness=").append(printCompleteness);}

		if(true || printUnique!=default_printUnique){sb.append(" printUMatches=").append(printUnique);}
		if(true || printUnique2!=default_printUnique2){sb.append(" printUnique2=").append(printUnique2);}
		if(true || printUnique3!=default_printUnique3){sb.append(" printUnique3=").append(printUnique3);}
		if(true || printUContam!=default_printUContam){sb.append(" printUContam=").append(printUContam);}
		if(true || printNoHit!=default_printNoHit){sb.append(" printNoHit=").append(printNoHit);}
		if(true || printContam!=default_printContam){sb.append(" contam=").append(printContam);}
		if(true){sb.append(" contam2=").append(printContam2 ? TaxTree.extendedToLevel(contamLevel)+"" : "f");}
		
		if(true || printScore!=default_printScore){sb.append(" printScore=").append(printScore);}
		
		if(true || printDepth!=default_printDepth){sb.append(" printDepth=").append(printDepth);}
		if(true || printDepth2!=default_printDepth2){sb.append(" printDepth2=").append(printDepth2);}
		if(true || printActualDepth!=default_printActualDepth){sb.append(" printActualDepth=").append(printActualDepth);}
		if(true || printVolume!=default_printVolume){sb.append(" printVolume=").append(printVolume);}
		
		if(true || printMatches!=default_printMatches){sb.append(" printMatches=").append(printMatches);}
		if(true || printLength!=default_printLength){sb.append(" printLength=").append(printLength);}
		if(true || printTaxID!=default_printTaxID){sb.append(" printTaxID=").append(printTaxID);}
		if(true || printGSize!=default_printGSize){sb.append(" printGSize=").append(printGSize);}
		if(true || printGKmers!=default_printGKmers){sb.append(" printGKmers=").append(printGKmers);}
		if(true || printTaxName!=default_printTaxName){sb.append(" printTaxName=").append(printTaxName);}
		if(true || printGSeqs!=default_printGSeqs){sb.append(" printGSeqs=").append(printGSeqs);}
		if(true || printGBases!=default_printGBases){sb.append(" printGBases=").append(printGBases);}
		if(true || minEntropy!=default_minEntropy){sb.append(" minEntropy=").append(String.format("%.4f", minEntropy));}
		if(comparator!=Comparison.scoreComparator){sb.append(" ").append(comparator.toString());}
		
		if(taxFilter!=null){
			sb.append(" taxfilterlevel=").append(taxFilter.taxLevel());
			sb.append(" taxfilterinclude=").append(taxFilter.include());
			if(taxFilter.containsString()!=null){
				sb.append(" taxfilterstring=").append(taxFilter.containsString());
			}
			Integer[] temp=taxFilter.taxSet();
			sb.append(" taxfilterset=");
			if(temp!=null){
				for(int i=0; i<temp.length; i++){
					if(i>0){sb.append(',');}
					sb.append(temp[i]);
				}
			}
		}
		
		if(useTaxidName){sb.append(" useTaxidName=").append(useTaxidName);}
		if(useImgName){sb.append(" useImgName=").append(useImgName);}
		if(useTaxName){sb.append(" useTaxName=").append(useTaxName);}
		
		if(true){sb.append(" colors=").append(printColors ? TaxTree.extendedToLevel(colorLevel)+"" : "f");}
		
		if(minKeyOccuranceCount!=default_minKeyOccuranceCount){sb.append(" minKeyOccuranceCount=").append(minKeyOccuranceCount);}
		
//		if(printColors && colorLevel!=default_colorLevel){sb.append(" colorLevel=").append(TaxTree.extendedToLevel(colorLevel));}
		

		if(printRefDivisor){sb.append(" printRefDivisor=").append(printRefDivisor);}
		if(printQueryDivisor){sb.append(" printQueryDivisor=").append(printQueryDivisor);}
		if(printRefSize){sb.append(" printRefSize=").append(printRefSize);}
		if(printQuerySize){sb.append(" printQuerySize=").append(printQuerySize);}
		if(printContamHits){sb.append(" printContamHits=").append(printContamHits);}
		if(printIntersection){sb.append(" printIntersection=").append(printIntersection);}
		
		if(reads>-1){sb.append(" reads=").append(reads);}
		if(mode!=default_mode){sb.append(" mode=").append(mode);}
		if(samplerate!=default_samplerate){sb.append(" samplerate=").append(String.format(Locale.ROOT, "%.4f",samplerate));}
		
		
		sb.append('\n');
		return sb.toString();
	}
	
	public boolean compatible(){
		return SketchObject.k==k && SketchObject.k2==k2 && SketchObject.amino==amino && hashVersion==SketchObject.HASH_VERSION;
	}
	
	public void setPrintAll(){
		printTax=true;
		printOriginalName=true;
		printImg=true;
		printAni=true;
		printCompleteness=true;
		printScore=true;
		printDepth=true;
		printDepth2=true;
		printVolume=true;
		
		printMatches=true;
		printLength=true;
		printTaxID=true;
		printGSize=true;
		printGKmers=true;
		printTaxName=true;
		printGSeqs=true;
		printGBases=true;
		
//		printColors=true;

		printUnique=true;
		printUnique2=true;
		printUnique3=true;
		printUContam=true;
		printNoHit=true;
		printContam=true;
		printContam2=true;
		
		printRefDivisor=true;
		printQueryDivisor=true;
		printRefSize=true;
		printQuerySize=true;
		printContamHits=true;
	}
	

	
	/*--------------------------------------------------------------*/
	/*----------------          Formatting          ----------------*/
	/*--------------------------------------------------------------*/

	StringBuilder queryHeader(Sketch sk){
		StringBuilder sb=new StringBuilder();
		if(format>2){return sb;}
		
		String color=toColor(sk.taxID);
		if(color!=null){sb.append(color);}
		
		sb.append("\nQuery: ").append(sk.name()==null ? "." : sk.name());
		if(dbName!=null){sb.append("\tDB: ").append(dbName);}
		sb.append("\tSketchLen: ").append(sk.length());
		sb.append("\tSeqs: ").append(sk.genomeSequences).append(' ');
		sb.append("\tBases: ").append(sk.genomeSizeBases);
		sb.append("\tgSize: ").append(sk.genomeSizeEstimate());
		if(sk.probCorrect<1 && sk.probCorrect>0){sb.append("\tQuality: ").append(String.format("%.4f", sk.probCorrect));}
		if(sk.counts!=null){
			double d=Tools.averageDouble(sk.counts);
			sb.append("\tAvgCount: ").append(String.format("%.3f", d));
			sb.append("\tDepth: ").append(String.format("%.3f", Tools.observedToActualCoverage(d)));
		}

		if(sk.imgID>0){sb.append("\tIMG: ").append(sk.imgID);}
		if(sk.spid>0){sb.append("\tspid: ").append(sk.spid);}
		if(sk.taxID>0 && sk.taxID<SketchObject.minFakeID){sb.append("\tTaxID: ").append(sk.taxID);}

		if(printFileName && sk.fname()!=null && !sk.fname().equals(sk.name())){sb.append("\tFile: "+sk.fname());}
		if(printOriginalName && sk.name0()!=null && !sk.name0().equals(sk.name())){sb.append("\tSeqName: "+sk.name0());}
		
		if(sk.meta!=null){
			for(String st : sk.meta){
				sb.append("\t").append(st.replaceFirst(":", ": "));
			}
		}
		
		if(color!=null){sb.append(Colors.RESET);}
		
		return sb;
	}
	
	int toColorTid(final int taxID){
		if(!printColors || SketchObject.taxtree==null || taxID<=0 || taxID>=SketchObject.minFakeID){return 0;}
		TaxNode tn=SketchObject.taxtree.getNode(taxID);
		while(tn!=null && tn.id!=tn.pid && tn.levelExtended<colorLevel){
			tn=SketchObject.taxtree.getNode(tn.pid);
//			System.err.println(tn);
		}
		return tn==null || tn.levelExtended>=TaxTree.LIFE_E || (tn.levelExtended>colorLevel && tn.levelExtended>TaxTree.PHYLUM_E) ? 0 : tn.id;
	}
	
	String toColor(final int taxID){
		if(!printColors || SketchObject.taxtree==null || taxID<=0 || taxID>=SketchObject.minFakeID){return null;}
		TaxNode tn=SketchObject.taxtree.getNode(taxID);
		while(tn!=null && tn.id!=tn.pid && tn.levelExtended<colorLevel){
			tn=SketchObject.taxtree.getNode(tn.pid);
//			System.err.println(tn);
		}
		if(tn==null){
			return null;
		}else{
			if(tn.levelExtended>=TaxTree.LIFE_E || (tn.levelExtended>colorLevel && tn.levelExtended>TaxTree.PHYLUM_E)){return Colors.WHITE;}
			else{
//				System.err.println("*"+tn.id+", "+tn.id%Colors.colorArray.length);
				return Colors.colorArray[tn.id%Colors.colorArray.length];
			}
		}
	}
	
	String header(){
		if(format==3){return "Query\tReference\tANI";}
		
		StringBuilder sb=new StringBuilder();
		
		//Numeric fields
		if(true){sb.append("WKID");}
		if(true){sb.append("\tKID");}
		if(printAni){sb.append("\tANI");}
		if(printCompleteness){sb.append("\tComplt");}
		if(printContam){sb.append("\tContam");}
		if(printContam2){sb.append("\tContam2");}
		if(printUContam){sb.append("\tuContam");}
		if(printScore){sb.append("\tScore");}
		if(printDepth){sb.append("\tDepth");}
		if(printDepth2){sb.append("\tDepth2");}
		if(printVolume){sb.append("\tVolume");}
		if(printMatches){sb.append("\tMatches");}
		if(printUnique){sb.append("\tUnique");}
		if(printUnique2){sb.append("\tUnique2");}
		if(printUnique3){sb.append("\tUnique3");}
		if(printNoHit){sb.append("\tnoHit");}
		if(printLength){sb.append("\tLength");}
		if(printTaxID){sb.append("\tTaxID");}
		if(printImg){sb.append("\tImgID    ");}
		if(printGBases){sb.append("\tgBases");}
		if(printGKmers){sb.append("\tgKmers");}
		if(printGSize){sb.append("\tgSize");}
		if(printGSeqs){sb.append("\tgSeqs");}
		
		
		//Raw fields
		if(printRefDivisor){sb.append("\trDiv");}
		if(printQueryDivisor){sb.append("\tqDiv");}
		if(printRefSize){sb.append("\trSize");}
		if(printQuerySize){sb.append("\tqSize");}
		if(printContamHits){sb.append("\tcHits");}
		
		//Text fields
		if(printTaxName){sb.append("\ttaxName");}
		if(printOriginalName){sb.append("\tseqName");}
		if(printTax && SketchObject.taxtree!=null){sb.append("\ttaxonomy");}
		
		return sb.toString();
	}
	
	void formatComparisonColumnwise(Comparison c, StringBuilder sb, int prevTid){
		final int tid=c.taxID;
		boolean reset=false;
		
		if(printColors){
			final int ctid=toColorTid(tid);
			final int prevCtid=toColorTid(prevTid);

			final int cnum=ctid%Colors.colorArray.length;
			final int prevCnum=prevCtid%Colors.colorArray.length;

			String color=toColor(tid);
			String underline=(printColors && cnum==prevCnum && ctid!=prevCtid && (ctid>1 && prevCtid>1) ? Colors.UNDERLINE : null);

			if(color!=null){sb.append(color);}
			if(underline!=null){sb.append(underline);}
			reset=(color!=null || underline!=null);
			
//			System.err.println((color==null ? "" : color)+(underline==null ? "" : underline)+
//					tid+", "+prevTid+";     \t"+ctid+", "+prevCtid+";     \t"+cnum+", "+prevCnum+"; \t"+((underline!=null)+"")+Colors.RESET);
//			System.err.println(color==null ? "null" : color.substring(1));
		}
		
		sb.append(String.format(Locale.ROOT, "%.2f%%\t%.2f%%", 100*c.idMinDivisor(), 100*c.idMaxDivisor()));
		
		if(printAni){sb.append(String.format(Locale.ROOT, "\t%.2f%%", 100*c.ani()));}
		if(printCompleteness){sb.append(String.format(Locale.ROOT, "\t%.2f%%", 100*c.completeness()));}
		if(printContam){sb.append(String.format(Locale.ROOT, "\t%.2f%%", 100*c.contamFraction()));}
		if(printContam2){sb.append(String.format(Locale.ROOT, "\t%.2f%%", 100*c.contam2Fraction()));}
		if(printUContam){sb.append(String.format(Locale.ROOT, "\t%.2f%%", 100*c.uContamFraction()));}
		if(printScore){sb.append('\t').append(c.scoreS());}
		
		if(printDepth){sb.append('\t').append(c.depthS(printActualDepth));}
		if(printDepth2){sb.append('\t').append(c.depth2S(printActualDepth));}
		if(printVolume){sb.append('\t').append(c.volumeS());}
		
		if(printMatches){sb.append('\t').append(c.hits);}
		if(printUnique){sb.append('\t').append(c.uHits());}
		if(printUnique2){sb.append('\t').append(c.unique2);}
		if(printUnique3){sb.append('\t').append(c.unique3);}
		if(printNoHit){sb.append('\t').append(c.noHits);}
		if(printLength){sb.append('\t').append( c.maxDivisor());}
		if(printTaxID){sb.append('\t').append(tid>=SketchObject.minFakeID ? -1 : tid);}
		if(printImg){sb.append(String.format(Locale.ROOT, "\t%d", c.imgID()));}
		if(printGBases){sb.append('\t').append(c.genomeSizeBases());}
		if(printGKmers){sb.append('\t').append(c.genomeSizeKmers());}
		if(printGSize){sb.append('\t').append(c.genomeSizeEstimate());}
		if(printGSeqs){sb.append('\t').append(c.genomeSequences());}
		
		//Raw fields
		if(printRefDivisor){sb.append('\t').append(c.refDivisor);}
		if(printQueryDivisor){sb.append('\t').append(c.queryDivisor);}
		if(printRefSize){sb.append('\t').append(c.refSize);}
		if(printQuerySize){sb.append('\t').append(c.querySize);}
		if(printContamHits){sb.append('\t').append(c.contamHits);}
		
		//Text fields
		if(printTaxName){sb.append('\t').append(c.taxName()==null ? "." : c.taxName());}
		if(printOriginalName || (c.taxName()==null && c.name0()!=null)){sb.append('\t').append(c.name0()==null ? "." : c.name0());}
		if(printTax && SketchObject.taxtree!=null){
			sb.append('\t');
			TaxNode tn=null;
			if(tid>0 && tid<SketchObject.minFakeID){
				tn=SketchObject.taxtree.getNode(tid);
			}

			if(tn!=null){
				sb.append(SketchObject.taxtree.toSemicolon(tn, SketchObject.skipNonCanonical));
			}else{
				sb.append('.');
			}
		}
		
		if(reset){sb.append(Colors.RESET);}
		
		sb.append('\n');
		
		if(printIntersection){
			Sketch intersection=Sketch.intersection(c.a, c.b);
			sb.append(intersection.toString());
			sb.append('\n');
		}
		
	}
	
	void formatComparison3Column(Comparison c, StringBuilder sb, int prevTid, Sketch query){

		final String qName=useTaxidName ? ""+query.taxID : useImgName ?  ""+query.imgID : useTaxName ? ""+query.taxName() : query.name();
		final String rName=useTaxidName ? ""+c.taxID() : useImgName ?  ""+c.imgID() : useTaxName ? ""+c.taxName() : c.name();
		final int tid=c.taxID;
		boolean reset=false;
		
		sb.append(qName).append('\t');
		
		if(printColors){
			final int ctid=toColorTid(tid);
			final int prevCtid=toColorTid(prevTid);

			final int cnum=ctid%Colors.colorArray.length;
			final int prevCnum=prevCtid%Colors.colorArray.length;

			String color=toColor(tid);
			String underline=(printColors && cnum==prevCnum && ctid!=prevCtid && (ctid>1 && prevCtid>1) ? Colors.UNDERLINE : null);

			if(color!=null){sb.append(color);}
			if(underline!=null){sb.append(underline);}
			reset=(color!=null || underline!=null);
			
//			System.err.println((color==null ? "" : color)+(underline==null ? "" : underline)+
//					tid+", "+prevTid+";     \t"+ctid+", "+prevCtid+";     \t"+cnum+", "+prevCnum+"; \t"+((underline!=null)+"")+Colors.RESET);
//			System.err.println(color==null ? "null" : color.substring(1));
		}
		
		sb.append(rName).append('\t').append(String.format(Locale.ROOT, "\t%.2f", 100*c.ani()));
		
		if(reset){sb.append(Colors.RESET);}
		
		sb.append('\n');
	}
	
	void formatComparison(Comparison c, StringBuilder sb, int prevTaxID, Sketch query){
		if(format==2){
			formatComparisonColumnwise(c, sb, prevTaxID);
			return;
		}else if(format==3){
			formatComparison3Column(c, sb, prevTaxID, query);
			return;
		}
		String complt=(printCompleteness ? String.format(Locale.ROOT, "\tcomplt %.2f%%%%", 100*c.completeness()) : "");
		String contam=(printContam ? String.format(Locale.ROOT, "\tcontam %.2f%%%%", 100*c.contamFraction()) : "");
//		String score=(printScore ? String.format(Locale.ROOT, "\tscore %.2f", c.score2()) : "");
		String score=(printScore ? "\tscore "+c.scoreS() : "");
		String depth=(printDepth ? "\tdepth "+c.depthS(printActualDepth) : "");
		String depth2=(printDepth2 ? "\tdepth2 "+c.depth2S(printActualDepth) : "");
		String volume=(printVolume ? "\tvolume "+c.volumeS() : "");
		String ccs=complt+contam+score;
		
		if(format==0){
			sb.append(String.format(Locale.ROOT, "WKID %.2f%%\tKID %.2f%%"+ccs+"\tmatches %d\tcompared %d",
					100*c.idMinDivisor(), 100*c.idMaxDivisor(), c.hits, c.minDivisor())+"\ttaxID "+c.taxID()+
					(printImg ? "\timgID "+c.imgID() : "")+"\tgKmers "+c.genomeSizeKmers()+"\t"+
					(c.taxName()==null ? "." : c.taxName())+
					((printOriginalName || (c.taxName()==null && c.name0()!=null)) ? "\t"+(c.name0()==null ? "." : c.name0()) : "")+"\n");
			if(printTax && SketchObject.taxtree!=null){
				if(c.taxID()>=0 && c.taxID()<SketchObject.minFakeID){
					TaxNode tn=SketchObject.taxtree.getNode(c.taxID());
					if(tn!=null){
						PrintTaxonomy.printTaxonomy(tn, sb, SketchObject.taxtree, TaxTree.DOMAIN, SketchObject.skipNonCanonical);
					}
				}
				sb.append('\n');
			}
		}else{
			ArrayList<TaxNode> tnl=new ArrayList<TaxNode>();
			if(SketchObject.taxtree!=null && c.taxID()>=0 && c.taxID()<SketchObject.minFakeID){
				TaxNode tn=SketchObject.taxtree.getNode(c.taxID());
				while(tn!=null && tn.pid!=tn.id && tn.level<=TaxTree.DOMAIN){
					tnl.add(tn);
					tn=SketchObject.taxtree.getNode(tn.pid);
				}
			}
			
			sb.append(String.format(Locale.ROOT, "WKID %.2f%%\tKID %.2f%%"+ccs+"\tmatches %d\tcompared %d\t",
					100*c.idMinDivisor(), 100*c.idMaxDivisor(), c.hits, c.minDivisor()));
			sb.append("\ttaxID ").append(c.taxID()).append('\t');
			if(printImg){sb.append("\timgID ").append(c.imgID()).append('\t');}
			sb.append(c.taxName()).append('\t');
			if(printOriginalName || (c.taxName()==null && c.name0()!=null)){sb.append(c.name0()).append('\t');}
			
			if(printTax){
				for(int i=tnl.size()-1; i>=0; i--){
					TaxNode tn=tnl.get(i);
					sb.append(tn.name);
					if(i>0){sb.append(';');}
				}
			}
			sb.append('\n');
			
			tnl.clear();
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	//These are shared with SketchObject
	//They do not affect anything and are just for the server to validate remote settings.
	private int hashVersion=SketchObject.HASH_VERSION;
	private int k=SketchObject.k;
	private int k2=SketchObject.k2;
	private boolean amino=SketchObject.amino;
	
	boolean postParsed=false;
	
	boolean amino(){return amino;}
	
	//These are unique
	public int maxRecords=default_maxRecords;
	public float minANI=0;
	public float minWKID=default_minWKID;
	public int format=default_format;
	
	public int minHits=default_minHits;
	public int taxLevel=default_taxLevel;
	public int mode=default_mode;
	public float samplerate=default_samplerate;
	public long reads=default_reads;
	public int minKeyOccuranceCount=default_minKeyOccuranceCount;
	
	public String dbName=null;
	
	/*--------------------------------------------------------------*/
	/*----------------         Print Columns        ----------------*/
	/*--------------------------------------------------------------*/
	
	//For format 2
	public boolean printTax=default_printTax;
	public boolean printOriginalName=default_printOriginalName;
	public boolean printFileName=default_printFileName;
	public boolean printImg=default_printImg;
	public boolean printAni=default_printAni;
	public boolean printCompleteness=default_printCompleteness;
	public boolean printScore=default_printScore;

	private boolean trackCounts=default_trackCounts;
	public boolean printDepth=default_printDepth;
	public boolean printDepth2=default_printDepth2;
	public boolean printActualDepth=default_printActualDepth;
	public boolean printVolume=default_printVolume;
	
	public boolean printLength=default_printLength;
	public boolean printTaxID=default_printTaxID;
	public boolean printGSize=default_printGSize;
	public boolean printGKmers=default_printGKmers;
	public boolean printTaxName=default_printTaxName;
	public boolean printGSeqs=default_printGSeqs;
	public boolean printGBases=default_printGBases;
	
	public float minEntropy=default_minEntropy;

	public boolean printUnique=default_printUnique;
	public boolean printUnique2=default_printUnique2;
	public boolean printUnique3=default_printUnique3;
	public boolean printUContam=default_printUContam;
	public boolean printNoHit=default_printNoHit;

	public boolean printColors=default_printColors;
	public boolean setColors=false;
	public int colorLevel=default_colorLevel;
	
	/** TODO: Note this is conflated between printing %contam and calculating things based on contam hits. */
	public boolean printContam=default_printContam;
	public boolean printContam2=default_printContam2;
	private int contamLevel=default_contamLevel;
	
	/** Raw fields */
	public boolean printMatches=default_printMatches;
	
	public boolean printRefDivisor=false;
	public boolean printQueryDivisor=false;
	public boolean printRefSize=false;
	public boolean printQuerySize=false;
	public boolean printContamHits=false;

	public boolean printIntersection=false;
	
	//For format 3
	public boolean useTaxidName=false;
	public boolean useImgName=false;
	public boolean useTaxName=false;
	
	public TaxFilter taxFilter=null;

	/** Make sure the settings are consistent, for CompareSketch.
	 * This is not yet complete. */
	public boolean checkValid(){
		if(printUnique2 || printUnique3){
			assert(contamLevel()>=TaxTree.SUBSPECIES_E);
			assert(needContamCounts());
			assert(SketchObject.makeIndex);
			assert(SketchObject.taxtree!=null);
		}
		if(printContam2){
			assert(contamLevel()>=TaxTree.SUBSPECIES_E);
			assert(needContamCounts());
			assert(SketchObject.makeIndex);
			assert(SketchObject.taxtree!=null);
		}
		return true;
	}
	
	public boolean trackCounts() {
		return trackCounts || printDepth || printDepth2 || printVolume || comparator!=Comparison.scoreComparator; //|| minKeyOccuranceCount>1;
	}
	
	public boolean needContamCounts() {
		return printContam || printContam2 || printContamHits || printUnique || printUnique2 || printUnique3 || printUContam || printNoHit; // || true
	}
	
	public boolean needIndex(){
		return printContam2 || printUnique2 || printUnique3;
	}

	public int contamLevel() {
		return needIndex() ? contamLevel : -1;
	}
	
	public int compare(Comparison a, Comparison b){
		return comparator.compare(a, b);
	}
	
	public Comparator<Comparison> comparator=Comparison.scoreComparator;
	
	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/
	
	public static final int default_maxRecords=20;
	public static final float default_minWKID=0.0001f;
	public static final int default_format=2;
	public static final boolean default_printTax=false;
	public static final boolean default_printOriginalName=false;
	public static final boolean default_printFileName=false;
	public static final boolean default_printImg=false;
	public static final boolean default_printAni=true;
	public static final boolean default_printCompleteness=true;
	public static final boolean default_printScore=false;
	
	public static final boolean default_trackCounts=false;
	public static final boolean default_printDepth=false;
	public static final boolean default_printDepth2=false;
	public static final boolean default_printActualDepth=true;
	public static final boolean default_printVolume=false;

	public static final boolean default_printContam=true;
	public static final boolean default_printContam2=false;
	
	public static final boolean default_printMatches=true;
	public static final boolean default_printLength=false;
	public static final boolean default_printTaxID=true;
	public static final boolean default_printGSize=true;
	public static final boolean default_printGKmers=false;
	public static final boolean default_printTaxName=true;
	public static final boolean default_printGSeqs=true;
	public static final boolean default_printGBases=false;

	public static final float default_minEntropy=0.66f;

	public static final boolean default_printUnique=true;
	public static final boolean default_printUnique2=false;
	public static final boolean default_printUnique3=false;
	public static final boolean default_printUContam=false;
	public static final boolean default_printNoHit=true;

	public static final boolean default_printColors=true;
	public static final int default_colorLevel=TaxTree.FAMILY_E;

	public static final int default_taxLevel=TaxTree.SPECIES;
	public static final int default_contamLevel=TaxTree.GENUS_E;
	
	public static final int default_mode=SketchObject.ONE_SKETCH;
	
	public static final int default_minHits=3;
	public static final float default_samplerate=1;
	public static final long default_reads=-1;
	public static final int default_minKeyOccuranceCount=1;
	
}
