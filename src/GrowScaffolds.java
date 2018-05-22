import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.lang.ProcessBuilder.Redirect;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Map.Entry;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class GrowScaffolds {
	private Map<String, Integer> contigToIndexMap;
	private Map<Integer, String> indexToContigMap;
	private Map<Integer, Boolean> contigIds;
	private Map<Integer, Double> contigAbundance;	
	private List<Set<Integer>> clusters;
	private Set<Integer> goodCovContigs;	
	private Map<Integer, Double> binAbundanceMap;
	private boolean ASC = true;
    private boolean DESC = false;    
	private String read1 = "";
	private String read2 = "";	
	private String outputDirName = "";
	private String FVEResDir = "";
	private String dbType = "";
	private String dbDir = "";
	private String dbFastaFile = "";
	private String dbFaiFile = "";
	private String dbClusterFile = "";
	private int minFoldCov = 3;
	private int minScaffoldLen = 2000;
	private int topBins = 100;
	private double avgReadLen = 0;
	private int numThreads = 4;
	private double cdHitCuttOff = 0.95;
	private String spadesKmerlen = "default";
	
	private void setDbFiles() {
		if (dbType.equals("gov")) {
			dbFastaFile = "gov_viral_populations.fasta";
			dbFaiFile = "gov_viral_populations.fasta.fai";
			dbClusterFile = "gov_viral_populations_clusters.txt";
		} 
		else if (dbType.equals("gov_epi_mes")) {
            dbFastaFile = "GOV_viral_contigs_EPI_MES.fna";
            dbFaiFile = "GOV_viral_contigs_EPI_MES.fna.fai";
            dbClusterFile = "GOV_viral_contigs_EPI_MES_clusters.txt";
        }
		else if (dbType.equals("imgvr")) {
			dbFastaFile = "mVCs_PaezEspino_Nature.fna";
			dbFaiFile = "mVCs_PaezEspino_Nature.fna.fai";
			dbClusterFile = "mVCs_PaezEspino_Nature_clusters.txt";
		}
	}
	
	private void generateCoverageStatForFVE() {
		String cmd = "";
		
		try {
			cmd = "cd " + FVEResDir + "\n"
					+ "bash pileup.sh "
					+ "in=FastViromeExplorer-reads-mapped-sorted.sam out=coverage-stat.txt twocolumn=t\n";
			
			FileWriter shellFileWriter = new FileWriter(outputDirName + "/run.sh");
			shellFileWriter.write("#!/bin/bash\n");
			shellFileWriter.write(cmd);
			shellFileWriter.close();

			ProcessBuilder builder = new ProcessBuilder("sh", outputDirName + "/run.sh");
			builder.redirectError(Redirect.appendTo(new File(outputDirName + "/log.txt")));
			Process process = builder.start();
			BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
			while (reader.readLine() != null) {
			}
			process.waitFor();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	private void getReadLen() {
		String cmd = "";
		try {          
            if (read1.endsWith(".gz")) {
            	cmd = "gzip -dc " + read1 
                        + " | awk 'NR%4 == 2 {lenSum+=length($0); readCount++;} END {print lenSum/readCount}'";
            } else {
            	cmd = "awk 'NR%4 == 2 {lenSum+=length($0); readCount++;} END {print lenSum/readCount}' "
                        + read1;
            }
            
            FileWriter shellFileWriter = new FileWriter(outputDirName + "/run.sh");
            shellFileWriter.write("#!/bin/bash\n");
            shellFileWriter.write(cmd);
            shellFileWriter.close();

            ProcessBuilder builder = new ProcessBuilder("sh", outputDirName + "/run.sh");
            builder.redirectError(Redirect.appendTo(new File(outputDirName + "/log.txt")));
            Process process = builder.start();
            BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
            String str = "";
            while ((str = reader.readLine()) != null) {
                avgReadLen = Double.parseDouble(str);
            }
            process.waitFor();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
	    
	    if (avgReadLen == 0) {
            System.out.println("Could not extract average read length from read file.");
            System.exit(1);
        }
	    System.out.println("Average estimated read length: " + avgReadLen);
	}
	
	private void createHashMap() {
		contigToIndexMap = new HashMap<String, Integer>();
		indexToContigMap = new HashMap<Integer, String>();
		BufferedReader br = null;
		try {
			br = new BufferedReader(new FileReader(dbDir + "/" + dbFaiFile));
			
			String str = "";
			String[] results;
			String contigId;
			int i = 0;
			while ((str = br.readLine()) != null) {
				str = str.trim();
				results = str.split("\t");
				contigId = results[0].trim();
				contigToIndexMap.put(contigId, i);
				indexToContigMap.put(i, contigId);
				i++;
			}
			br.close();
			//System.out.println("Map size: " + contigToIndexMap.size());
			//System.out.println("Map size: " + indexToContigMap.size());
		}
		catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	private void readFVERes() {
		contigIds = new HashMap<Integer, Boolean>();
		contigAbundance = new HashMap<Integer, Double>();
		BufferedReader br = null;
		try {
			br = new BufferedReader(new FileReader(FVEResDir 
					+ "/FastViromeExplorer-final-sorted-abundance.tsv"));
			
			String str = "";
			String[] results;
			String contigId;
			int index;
			double abundance;
			
			br.readLine();
			while ((str = br.readLine()) != null) {
				str = str.trim();
				results = str.split("\t");
				contigId = results[0].trim();
				abundance = Double.parseDouble(results[3].trim());
				index = contigToIndexMap.get(contigId);
				contigIds.put(index, false);
				contigAbundance.put(index, abundance);
			}
			br.close();
			//System.out.println("Num of all contigs: " + contigIds.size());
			//System.out.println("Num of all contigs: " + contigAbundance.size());
		}
		catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	private void getClusters() {
		clusters = new ArrayList<Set<Integer>>();
		BufferedReader br = null;
		try {
			br = new BufferedReader(new FileReader(dbDir + "/" + dbClusterFile));
			
			String str = "";
			String[] results;
			Set<Integer> cluster;
			while ((str = br.readLine()) != null) {
				str = str.trim();
				results = str.split("\t");
				cluster = new HashSet<Integer>();
				for (int i = 0; i < results.length; i++) {
					int index = Integer.parseInt(results[i]);
					if(contigIds.containsKey(index)) {
						cluster.add(index);
						if (contigIds.get(index)) {
							System.out.println("Something is wrong for: " + index);
						}
						contigIds.put(index, true);
					}
				}
				if(cluster.size() > 0) {
					clusters.add(cluster);
				}
			}
			br.close();
			for (int key:contigIds.keySet()) {
				if(!contigIds.get(key)) {
					cluster = new HashSet<Integer>();
					cluster.add(key);
					clusters.add(cluster);
				}
			}			
			contigIds = null;
		}
		catch (Exception e) {
			e.printStackTrace();
		}
		System.out.println("Total num of initial bins: " + clusters.size());
	}
	
	private boolean isCovGood(Set<Integer> cluster) {
		for (int i:cluster) {
			if(goodCovContigs.contains(i))
				return true;
		}
		return false;
	}
	
	private void checkCoverage() {
		goodCovContigs = new HashSet<Integer>();
		BufferedReader br = null;
		try {
			br = new BufferedReader(new FileReader(FVEResDir + "/coverage-stat.txt"));
			
			String str = "";
			String[] results;
			br.readLine();
			while ((str = br.readLine()) != null) {
				str = str.trim();
				results = str.split("\t");
				if(Double.parseDouble(results[1]) >= minFoldCov) {
					goodCovContigs.add(contigToIndexMap.get(results[0]));
				}
			}
			br.close();			
		}
		catch (Exception e) {
			e.printStackTrace();
		}
		
		//System.out.println("Good coverage contigs: " + goodCovContigs.size());
		
		List<Set<Integer>> newClusters = new ArrayList<Set<Integer>>();
		for (int i = 0; i < clusters.size(); i++) {
			Set<Integer> cluster = clusters.get(i);
			if (isCovGood(cluster)) {
				newClusters.add(cluster);
			}
		}
		clusters = newClusters;
		System.out.println("Total num of bins after checking coverage: " + clusters.size());

		goodCovContigs = null;
	}
	
	// sort a map
    private Map<Integer, Double> sortByComparator(Map<Integer, Double> unsortMap, final boolean order) {
        List<Entry<Integer, Double>> list = new LinkedList<Entry<Integer, Double>>(unsortMap.entrySet());

        // Sorting the list based on values
        Collections.sort(list, new Comparator<Entry<Integer, Double>>() {
            public int compare(Entry<Integer, Double> o1, Entry<Integer, Double> o2) {
                if (order) {
                    return o1.getValue().compareTo(o2.getValue());
                } else {
                    return o2.getValue().compareTo(o1.getValue());

                }
            }
        });

        // Maintaining insertion order with the help of LinkedList
        Map<Integer, Double> sortedMap = new LinkedHashMap<Integer, Double>();
        int index = 0;
        for (Entry<Integer, Double> entry : list) {
            sortedMap.put(entry.getKey(), entry.getValue());
            index++;
            if(index == topBins)
            	break;
        }

        return sortedMap;
    }
	
	private void getTopBins() {
		binAbundanceMap = new HashMap<Integer, Double>();
		for (int i = 0; i < clusters.size(); i++) {
			Set<Integer> cluster = clusters.get(i);
			if (cluster.size() > 0) {
				double binAbundance = 0;
				for (int j:cluster) {
					binAbundance += contigAbundance.get(j);
				}
				binAbundanceMap.put(i, binAbundance);
			}
		}
		binAbundanceMap = sortByComparator(binAbundanceMap, DESC);
		//System.out.println("Sorted map size: " + binAbundanceMap.size());
		List<Set<Integer>> newClusters = new ArrayList<Set<Integer>>();
		for (int key:binAbundanceMap.keySet()) {
			Set<Integer> cluster = clusters.get(key);
			newClusters.add(cluster);
		}

		clusters = newClusters;
		contigAbundance = null;
		binAbundanceMap = null;
		System.out.println("Total num of final bins after getting top bins: " + clusters.size());
		/*int elementsInClusters = 0;
		
		for (int i = 0; i < clusters.size(); i++) {
			Set<Integer> cluster = clusters.get(i);
			elementsInClusters += cluster.size();
		}
		System.out.println("Num of elements in final clusters: " + elementsInClusters);*/
	}
	
	public void writeFastqFilesAndAssembly() {		
		//call metaSpades in parallel
		ExecutorService executor = Executors.newFixedThreadPool(numThreads);
		for (int i = 0; i < clusters.size(); i++) {
			Runnable worker = new MyRunnableGenerateSeeds(i);
			executor.execute(worker);
		}
		executor.shutdown();
		//wait until all threads are finished
		while (!executor.isTerminated()) {
			
		}
	}
	
	public class MyRunnableGenerateSeeds implements Runnable {
		private final int binNum;
		
		public MyRunnableGenerateSeeds(int binNum) {
			this.binNum = binNum;
		}
		
		@Override
		public void run() {
			Set<Integer> cluster = clusters.get(binNum);
			if (cluster.size() > 0)
			{				
				BufferedReader br = null;
				BufferedWriter bw = null;
				try {
					br = new BufferedReader(new FileReader(FVEResDir + "/FastViromeExplorer-reads-mapped-sorted.sam"));
					
					File dir = new File(outputDirName + "/bin-" + binNum);
					dir.mkdir();
					bw = new BufferedWriter(new FileWriter(outputDirName + "/bin-" + binNum + "/readIds.txt"));
					
					String str = "";
					String[] results;
					int contigIndex = 0;
					while ((str = br.readLine()) != null) {
						if (!str.startsWith("@")) {
							str = str.trim();
							results = str.split("\t");
							contigIndex = contigToIndexMap.get(results[2]);
							if (cluster.contains(contigIndex)) {
								bw.write(results[0] + "\n");
							}
						}				
					}
					br.close();
					bw.close();
				}
				catch (Exception e) {
					e.printStackTrace();
				}
				
				String cmd = "";
				try 
				{
					if (read2.isEmpty()) {
					    cmd = "cd " + outputDirName + "\n"
                            + "cd bin-" + binNum + "\n"
                            + "bash filterbyname.sh "
                            + "in=" + read1 + " "
                            + "out=reads_1.fastq names=readIds.txt "
                            + "include=t\n"
                            + "spades.py "
                            + "-o metaSpades-res -s reads_1.fastq --only-assembler\n"
                            + "cd metaSpades-res\n"
                            + "samtools faidx scaffolds.fasta\n"
                            + "cd-hit-est -i scaffolds.fasta -o cd-hit-scaffolds.fasta "
                            + "-c " + cdHitCuttOff + " -n 10 -d 0 -M 16000 -T 8\n"
                            + "samtools faidx cd-hit-scaffolds.fasta\n";
					}
					else {
        				    cmd = "cd " + outputDirName + "\n"
        							+ "cd bin-" + binNum + "\n"
        							+ "bash filterbyname.sh "
        							+ "in1=" + read1 + " "
        							+ "in2=" + read2 + " "
        							+ "out1=reads_1.fastq out2=reads_2.fastq names=readIds.txt "
        							+ "include=t\n"
        							+ "spades.py --meta "
        							+ "-o metaSpades-res -1 reads_1.fastq -2 reads_2.fastq --only-assembler\n"
        							+ "cd metaSpades-res\n"
        							+ "samtools faidx scaffolds.fasta\n"
        							+ "cd-hit-est -i scaffolds.fasta -o cd-hit-scaffolds.fasta "
        							+ "-c " + cdHitCuttOff + " -n 10 -d 0 -M 16000 -T 8\n"
        							+ "samtools faidx cd-hit-scaffolds.fasta\n";
					}
					
					FileWriter shellFileWriter = new FileWriter(outputDirName + "/bin-" + binNum + "/run.sh");
					shellFileWriter.write("#!/bin/bash\n");
					shellFileWriter.write(cmd);
					shellFileWriter.close();

					ProcessBuilder builder = new ProcessBuilder("sh", outputDirName + "/bin-" + binNum + "/run.sh");
					builder.redirectError(Redirect.appendTo(new File(outputDirName + "/log.txt")));
					Process process = builder.start();
					BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
					while (reader.readLine() != null) {
					}
					process.waitFor();
				} catch (Exception ex) {
					ex.printStackTrace();
				}
				//System.out.println("Finished processing bin " + binNum);
			}
		}
	}
	
	public int getSeedScaffolds() {
		int totalBin = clusters.size();
		if (totalBin == 0) {
			return 0;
		}
		File dir = new File(outputDirName + "/seed-scaffolds");
		dir.mkdir();
		int scaffoldNum = 0;		
		BufferedReader br = null;
		BufferedWriter bw = null;
		int binsFailedAssembly = 0;
		
		for (int i = 0; i < totalBin; i++) {
			File file = new File(outputDirName + "/bin-" + i + "/metaSpades-res/cd-hit-scaffolds.fasta.fai");
			if (file.exists())
			{
				try 
				{
					br = new BufferedReader(new FileReader(file));
					bw = new BufferedWriter(new FileWriter(outputDirName + "/seed-scaffolds/binInfo.txt", true));
					
					String str = "";
					String[] results;
					int length = 0;
					while ((str = br.readLine()) != null) 
					{
						results = str.split("\t");
						length = Integer.parseInt(results[1].trim());
						if (length >= minScaffoldLen) 
						{
							String cmd = "cd " + outputDirName + "\n" 
									+ "cd seed-scaffolds\n"
									+ "mkdir scaffold-" + scaffoldNum + "\n"
									+ "bash filterbyname.sh "
									+ "in=../bin-" + i + "/metaSpades-res/cd-hit-scaffolds.fasta "
									+ "out=scaffold-" + scaffoldNum + "/scaffold.fasta names=" + results[0]
									+ " include=t\n"
									+ "cd scaffold-" + scaffoldNum + "\n"
									+ "samtools faidx scaffold.fasta\n";
							
							FileWriter shellFileWriter = new FileWriter(outputDirName + "/run.sh");
							shellFileWriter.write("#!/bin/bash\n");
							shellFileWriter.write(cmd);
							shellFileWriter.close();

							ProcessBuilder builder = new ProcessBuilder("sh", outputDirName + "/run.sh");
							builder.redirectError(Redirect.appendTo(new File(outputDirName + "/log.txt")));
							Process process = builder.start();
							BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
							while (reader.readLine() != null) {
							}
							process.waitFor();
							bw.write("scaffold-" + scaffoldNum + "\tLength: " + length + "\tbin-" + i + "\t" + clusters.get(i) + "\n");
							scaffoldNum++;
						}						
					}
					br.close();
					bw.close();
				}
				catch (Exception e) {
					e.printStackTrace();
				}
			}
			else {
				//System.out.println("No assembly for bin: " + i);
				binsFailedAssembly++;
			}
		}
		//System.out.println("Total bins failed in metaSpades assembly: " + binsFailedAssembly);
		//System.out.println("Total seed scaffolds: " + scaffoldNum);
		return scaffoldNum;
	}
	
	private void createBed(String outputDir) {
		//move contig.fasta and contig.fasta.fai from spades-res folder to outputDir
		File scaffoldFile = new File(outputDir + "/tmp/spades-res/scaffold.fasta");
		if (scaffoldFile.exists() && !scaffoldFile.isDirectory()) 
		{
			String cmd = "";
			try {
				cmd = "cd " + outputDir + "\n"
						+ "rm scaffold.fasta*\n"
						+ "cd tmp/spades-res\n" 
						+ "mv scaffold.fasta " + outputDir + "\n"
						+ "mv scaffold.fasta.fai " + outputDir + "\n";

				FileWriter shellFileWriter = new FileWriter(outputDir + "/run.sh");
				shellFileWriter.write("#!/bin/bash\n");
				shellFileWriter.write(cmd);
				shellFileWriter.close();

				ProcessBuilder builder = new ProcessBuilder("sh", outputDir + "/run.sh");
				builder.redirectError(Redirect.appendTo(new File(outputDir + "/log.txt")));
				Process process = builder.start();
				BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
				while (reader.readLine() != null) {
				}
				process.waitFor();
			} catch (Exception ex) {
				ex.printStackTrace();
			}			
		}		
		
		BufferedReader br = null;
		BufferedWriter bw = null;
		BufferedWriter bwOutLog = null;
		try {
			br = new BufferedReader(new FileReader(outputDir + "/scaffold.fasta.fai"));
			bw = new BufferedWriter(new FileWriter(outputDir + "/scaffold-start-end.bed"));
			bwOutLog = new BufferedWriter(new FileWriter(outputDir + "/output-log.txt", true));
			String str = br.readLine();
			String[] results = str.split("\t");
			int scaffoldLength = Integer.parseInt(results[1].trim());
			if (scaffoldLength >= avgReadLen * 2) {
				String scaffoldId = results[0].trim();
				bw.write(scaffoldId + "\t" + 0 + "\t" + (int) Math.ceil(avgReadLen * 1.5) + "\n");
				bw.write(scaffoldId + "\t" + (int) Math.ceil(scaffoldLength - avgReadLen * 1.5) + "\t" + scaffoldLength
						+ "\n");
				bwOutLog.write("Trying to grow scaffold " + scaffoldId + " with length "
						+ scaffoldLength + "\n");
			}
			br.close();
			bw.close();
			bwOutLog.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	private void runAlignment(String outputDir) {
		String cmd = "";
		try {
		    if (read2.isEmpty()) {
		        cmd = "cd " + outputDir + "\n"
                    + "bedtools getfasta -fi scaffold.fasta -bed scaffold-start-end.bed -fo scaffold-start-end.fasta\n"
                    + "rm -r salmon-index\n"
                    + "rm -r salmon-res\n"
                    + "rm salmon-mapped.sam\n"
                    + "salmon index -t scaffold-start-end.fasta -i salmon-index\n"
                    + "/usr/bin/time -f \"\t%E Elasped Real Time\" salmon quant -i salmon-index -l A "
                    + "-r " + read1 + " -o salmon-res --writeMappings -p 16 "
                    + "| samtools view -bS - | samtools view -h -F 0x04 - > salmon-mapped.sam\n";
		    }
		    else {
        			cmd = "cd " + outputDir + "\n"
        					+ "bedtools getfasta -fi scaffold.fasta -bed scaffold-start-end.bed -fo scaffold-start-end.fasta\n"
        					+ "rm -r salmon-index\n"
        					+ "rm -r salmon-res\n"
        					+ "rm salmon-mapped.sam\n"
        					+ "salmon index -t scaffold-start-end.fasta -i salmon-index\n"
        					+ "/usr/bin/time -f \"\t%E Elasped Real Time\" salmon quant -i salmon-index -l A "
        					+ "-1 " + read1 + " -2 " + read2 + " -o salmon-res --writeMappings -p 16 "
        					+ "| samtools view -bS - | samtools view -h -F 0x04 - > salmon-mapped.sam\n";
		    }
			
			FileWriter shellFileWriter = new FileWriter(outputDir + "/run.sh");
			shellFileWriter.write("#!/bin/bash\n");
			shellFileWriter.write(cmd);
			shellFileWriter.close();

			ProcessBuilder builder = new ProcessBuilder("sh", outputDir + "/run.sh");
			builder.redirectError(Redirect.appendTo(new File(outputDir + "/log-alignment.txt")));
			Process process = builder.start();
			BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
			while (reader.readLine() != null) {
			}
			process.waitFor();
		} catch (Exception ex) {
			ex.printStackTrace();
		}
	}
	
	private void getMappedReads(String outputDir) {		
		String cmd = "";
		try {
		    if (read2.isEmpty()) {
		        cmd = "cd " + outputDir + "\n" 
                    + "rm -r tmp\n" 
                    + "mkdir tmp\n"
                    + "bash filterbyname.sh in=" + read1
                    + " out=tmp/mapped_reads_1.fastq names="
                    + "salmon-mapped.sam include=t\n";
		    }
		    else {
		        cmd = "cd " + outputDir + "\n" 
                    + "rm -r tmp\n" 
                    + "mkdir tmp\n"
                    + "bash filterbyname.sh in=" + read1 + " in2=" + read2
                    + " out=tmp/mapped_reads_1.fastq out2=tmp/mapped_reads_2.fastq names="
                    + "salmon-mapped.sam include=t\n";
		    }			

			FileWriter shellFileWriter = new FileWriter(outputDir + "/run.sh");
			shellFileWriter.write("#!/bin/bash\n");
			shellFileWriter.write(cmd);
			shellFileWriter.close();

			ProcessBuilder builder = new ProcessBuilder("sh", outputDir + "/run.sh");
			builder.redirectError(Redirect.appendTo(new File(outputDir + "/log.txt")));
			Process process = builder.start();
			BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
			while (reader.readLine() != null) {
			}
			process.waitFor();
		} catch (Exception ex) {
			ex.printStackTrace();
		}
	}
	
	private void runSpades(String outputDir) {
		String cmd = "";
		try {
		    if (read2.isEmpty()) {
		        if (spadesKmerlen.equals("default")) {
		            cmd = "cd " + outputDir + "/tmp\n" 
	                    + "/usr/bin/time -f \"\t%E Elasped Real Time\" spades.py -o "
	                    + "spades-res -s mapped_reads_1.fastq "
	                    + "--trusted-contigs ../scaffold.fasta"
	                    + " --only-assembler\n"
	                    + "cd spades-res\n"
	                    + "samtools faidx scaffolds.fasta\n";
		        }
		        else {
        		        cmd = "cd " + outputDir + "/tmp\n" 
                            + "/usr/bin/time -f \"\t%E Elasped Real Time\" spades.py -o "
                            + "spades-res -s mapped_reads_1.fastq "
                            + "--trusted-contigs ../scaffold.fasta -k " + spadesKmerlen
                            + " --only-assembler\n"
                            + "cd spades-res\n"
                            + "samtools faidx scaffolds.fasta\n";
		        }
		    }
		    else {
		        if (spadesKmerlen.equals("default")) {
		            cmd = "cd " + outputDir + "/tmp\n" 
                        + "/usr/bin/time -f \"\t%E Elasped Real Time\" spades.py -o "
                        + "spades-res -1 mapped_reads_1.fastq -2 mapped_reads_2.fastq "
                        + "--trusted-contigs ../scaffold.fasta"
                        + " --only-assembler\n"
                        + "cd spades-res\n"
                        + "samtools faidx scaffolds.fasta\n";
		        }
		        else {
            			cmd = "cd " + outputDir + "/tmp\n" 
            					+ "/usr/bin/time -f \"\t%E Elasped Real Time\" spades.py -o "
            					+ "spades-res -1 mapped_reads_1.fastq -2 mapped_reads_2.fastq "
            					+ "--trusted-contigs ../scaffold.fasta -k " + spadesKmerlen
            					+ " --only-assembler\n"
            					+ "cd spades-res\n"
            					+ "samtools faidx scaffolds.fasta\n";
		        }
		    }

			FileWriter shellFileWriter = new FileWriter(outputDir + "/run.sh");
			shellFileWriter.write("#!/bin/bash\n");
			shellFileWriter.write(cmd);
			shellFileWriter.close();

			ProcessBuilder builder = new ProcessBuilder("sh", outputDir + "/run.sh");
			builder.redirectError(Redirect.appendTo(new File(outputDir + "/log-assembly.txt")));
			Process process = builder.start();
			BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
			while (reader.readLine() != null) {
			}
			process.waitFor();
		} catch (Exception ex) {
			ex.printStackTrace();
		}
	}
	
	private boolean firstExtend(String outputDir) {
		createBed(outputDir);
		runAlignment(outputDir);
		getMappedReads(outputDir);
		runSpades(outputDir);
		File scaffoldFile = new File(outputDir + "/tmp/spades-res/scaffolds.fasta");
		if (scaffoldFile.exists() && !scaffoldFile.isDirectory()) {
			return true;			
		}
		return false;
	}
	
	private int getScaffoldFromScaffolds(String outputDir) {	
		int maxScaffoldLength = 0;
		String maxScaffoldId = "";
		
		BufferedReader br = null;
		try {
			br = new BufferedReader(new FileReader(outputDir + "/tmp/spades-res/scaffolds.fasta.fai"));
			String str = br.readLine();
			str = str.trim();
			String[] results = str.split("\t");
			maxScaffoldId = results[0].trim();
			maxScaffoldLength = Integer.parseInt(results[1].trim());			
			br.close();
			
			if (!maxScaffoldId.equals("") && maxScaffoldLength != 0)
			{
				String cmd = "";
				
				cmd = "cd " + outputDir + "/tmp/spades-res\n"
						+ "bash filterbyname.sh "
						+ "in=scaffolds.fasta "
						+ "out=scaffold.fasta names=" + maxScaffoldId
						+ " include=t\n"
						+ "samtools faidx scaffold.fasta\n";
				
				FileWriter shellFileWriter = new FileWriter(outputDir + "/run.sh");
				shellFileWriter.write("#!/bin/bash\n");
				shellFileWriter.write(cmd);
				shellFileWriter.close();
	
				ProcessBuilder builder = new ProcessBuilder("sh", outputDir + "/run.sh");
				builder.redirectError(Redirect.appendTo(new File(outputDir + "/log.txt")));
				Process process = builder.start();
				BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
				while (reader.readLine() != null) {
				}
				process.waitFor();
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		return maxScaffoldLength;
	}
	
	private void growScaffold(int scaffoldNum) {
		int iteration = 1;
		String outputDir = outputDirName + "/seed-scaffolds/scaffold-" + scaffoldNum;
		boolean extendContig = firstExtend(outputDir);
		iteration++;
		int prevLength = minScaffoldLen;
		
		while (extendContig) {
			int currentLength = getScaffoldFromScaffolds(outputDir);
			if (currentLength > prevLength)
			{
				prevLength = currentLength;
				createBed(outputDir);
				runAlignment(outputDir);
				getMappedReads(outputDir);
				runSpades(outputDir);
				File scaffoldFile = new File(outputDir + "/tmp/spades-res/scaffolds.fasta");
				if (scaffoldFile.exists() && !scaffoldFile.isDirectory()) {
					iteration++;
					extendContig = true;
				} else {
					extendContig = false;
				}
			}
			else {
				extendContig = false;
			}
			
			if (extendContig && iteration > 1000) {
				extendContig = false;
			}
		}
	}
	
	public class MyRunnableGrowScaffold implements Runnable {
		private final int seedScaffoldNum;
		
		public MyRunnableGrowScaffold(int seedScaffoldNum) {
			this.seedScaffoldNum = seedScaffoldNum;
		}
		
		@Override
		public void run() {
			growScaffold(seedScaffoldNum);
		}
	}
	
	public void growSeedScaffolds(int totalSeedScaffolds, int seedExtensionStart) {		
		//extend contigs in parallel
		ExecutorService executor = Executors.newFixedThreadPool(numThreads);
		if (seedExtensionStart >= totalSeedScaffolds) {
		    System.out.println("Invalid seedExtension starting number. Changing it to 0");
		    seedExtensionStart = 0;
		}
		for (int i = seedExtensionStart; i < totalSeedScaffolds; i++) {
			Runnable worker = new MyRunnableGrowScaffold(i);
			executor.execute(worker);
		}
		executor.shutdown();
		//wait until all threads are finished
		while (!executor.isTerminated()) {
			
		}
		//System.out.println("Finished all threads");
	}
	
	private void getAllScaffolds(int totalSeedScaffolds) {		
		//copy all scaffold fasta to all-scaffolds.fasta
		File dir = new File(outputDirName + "/final-results");
		dir.mkdir();
		
		for (int i = 0; i < totalSeedScaffolds; i++) 
		{
			String scaffoldFileName = outputDirName
					+ "/seed-scaffolds/scaffold-" + i + "/scaffold.fasta";
			File scaffoldFile = new File(scaffoldFileName);
			if (scaffoldFile.exists() && !scaffoldFile.isDirectory()) 
			{
				String cmd = "";
				
				try {
					cmd = "cd " + outputDirName + "\n"
							+ "cat seed-scaffolds/scaffold-" + i 
							+ "/scaffold.fasta >> final-results/all-scaffolds.fasta\n";
					
					FileWriter shellFileWriter = new FileWriter(outputDirName + "/run.sh");
					shellFileWriter.write("#!/bin/bash\n");
					shellFileWriter.write(cmd);
					shellFileWriter.close();
	
					ProcessBuilder builder = new ProcessBuilder("sh", outputDirName + "/run.sh");
					builder.redirectError(Redirect.appendTo(new File(outputDirName + "/log.txt")));
					Process process = builder.start();
					BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
					while (reader.readLine() != null) {
					}
					process.waitFor();
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		}
		
		BufferedReader br = null;
		BufferedWriter bw = null;
		try {
			br = new BufferedReader(new FileReader(outputDirName + "/final-results/all-scaffolds.fasta"));
			bw = new BufferedWriter(new FileWriter(outputDirName + "/final-results/all-scaffolds-corrected.fasta"));
			String str = "";
			int scaffoldNum = 0;
			while ((str = br.readLine()) != null) 
			{
				str = str.trim();
				if (str.startsWith(">")) {
					str = ">" + scaffoldNum + "-" + str.substring(1);
					bw.write(str + "\n");
					scaffoldNum++;
				}
				else {
					bw.write(str + "\n");
				}
			}
			br.close();
			bw.close();
			//System.out.println("All scaffolds: " + scaffoldNum);
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		String cmd = "";
		
		try {
			cmd = "cd " + outputDirName + "/final-results\n"
					+ "bash sortbyname.sh "
					+ "in=all-scaffolds-corrected.fasta out=all-scaffolds-sorted.fasta length descending\n"
					+ "samtools faidx all-scaffolds-sorted.fasta\n";
			
			FileWriter shellFileWriter = new FileWriter(outputDirName + "/run.sh");
			shellFileWriter.write("#!/bin/bash\n");
			shellFileWriter.write(cmd);
			shellFileWriter.close();

			ProcessBuilder builder = new ProcessBuilder("sh", outputDirName + "/run.sh");
			builder.redirectError(Redirect.appendTo(new File(outputDirName + "/log.txt")));
			Process process = builder.start();
			BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
			while (reader.readLine() != null) {
			}
			process.waitFor();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	private void runCdhit() {
		String cmd = "";
		
		try {
			cmd = "cd " + outputDirName + "/final-results\n"
					+ "cd-hit-est -i all-scaffolds-sorted.fasta -o final-scaffolds.fasta -c " + cdHitCuttOff 
					+ " -n 10 -d 0 -M 16000 -T 8\n"
					+ "samtools faidx final-scaffolds.fasta\n";
			
			FileWriter shellFileWriter = new FileWriter(outputDirName + "/run.sh");
			shellFileWriter.write("#!/bin/bash\n");
			shellFileWriter.write(cmd);
			shellFileWriter.close();

			ProcessBuilder builder = new ProcessBuilder("sh", outputDirName + "/run.sh");
			builder.redirectError(Redirect.appendTo(new File(outputDirName + "/log.txt")));
			Process process = builder.start();
			BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
			while (reader.readLine() != null) {
			}
			process.waitFor();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public void getFinalScaffolds(int totalSeedScaffolds) {
		getAllScaffolds(totalSeedScaffolds);
		runCdhit();
	}
	
	private void getBinAndANIInfo() {
        Map<Integer, Integer> scaffoldToBinMap = new HashMap<Integer, Integer>();
        //get scaffold-bin map
        BufferedReader br = null;
        try {
            br = new BufferedReader(new FileReader(outputDirName + "/seed-scaffolds/binInfo.txt"));
            String str = "";
            String[] results;
            while ((str = br.readLine()) != null) 
            {
                str = str.trim();
                results = str.split("\t");
                int scaffold = Integer.parseInt(results[0].substring(results[0].indexOf('-')+1));
                int bin = Integer.parseInt(results[2].substring(results[2].indexOf('-')+1));
                scaffoldToBinMap.put(scaffold, bin);
            }
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        
        int originalScaffoldNum = 0;
        int binNum = 0;
        String fastaHeader = "";
        int sortedScaffoldNum = 0;
        Set<Integer> genomesInBin;
        String cmd = "";
        String ANIRes = "";
        double refANI = 0;
        double qryANI = 0;
        double alignedNuclRef = 0;
        double alignedNuclQry = 0;
        int refLength = 0;
        int qryLength = 0;
        
        BufferedReader brDnaDiff = null;
        BufferedWriter bwBin = null;
        BufferedWriter bwANI = null;
        try {
            br = new BufferedReader(new FileReader(outputDirName 
                    + "/final-results/final-scaffolds.fasta.fai"));
            bwBin = new BufferedWriter(new FileWriter(outputDirName 
                    + "/final-results/bin-info.tsv"));
            bwBin.write("ScaffoldNum\tScaffold\tBin(Genomes/contigs in this bin)\n");
            
            bwANI = new BufferedWriter(new FileWriter(outputDirName 
                    + "/final-results/ANI-info.tsv"));
            bwANI.write("ScaffoldNum\tScaffold\tdnadiff output from MUMmer tool\n");
            
            String str = "";
            String[] results;
            
            while ((str = br.readLine()) != null) 
            {
                str = str.trim();
                results = str.split("\t");
                fastaHeader = results[0].trim();
                
                originalScaffoldNum = Integer.parseInt(fastaHeader.substring(0, fastaHeader.indexOf('-')));
                binNum = scaffoldToBinMap.get(originalScaffoldNum);
                
                bwBin.write("Scaffold-" + sortedScaffoldNum + "\t" + fastaHeader + "\tBin-" + binNum + "(");
                genomesInBin = clusters.get(binNum);
                boolean firstElem = true;
                for (int genomeId:genomesInBin) 
                {
                    if (firstElem) {
                        bwBin.write(indexToContigMap.get(genomeId));
                        firstElem = false;
                    }
                    else {
                        bwBin.write(", " + indexToContigMap.get(genomeId));
                    }
                }
                bwBin.write(")\n");
                
                cmd = "";
                try {
                    cmd = "cd " + outputDirName + "/final-results\n"
                            + "mkdir scaffold-" + sortedScaffoldNum + "\n"
                            + "cd scaffold-" + sortedScaffoldNum + "\n"
                            + "bash filterbyname.sh in=../final-scaffolds.fasta "
                            + "out=scaffold-" + sortedScaffoldNum + ".fasta names=" + fastaHeader + " include=t\n"
                            + "samtools faidx scaffold-" + sortedScaffoldNum + ".fasta\n";
                    
                    FileWriter shellFileWriter = new FileWriter(outputDirName + "/final-results/run.sh");
                    shellFileWriter.write("#!/bin/bash\n");
                    shellFileWriter.write(cmd);
                    shellFileWriter.close();
    
                    ProcessBuilder builder = new ProcessBuilder("sh", outputDirName + "/final-results/run.sh");
                    builder.redirectError(Redirect.appendTo(new File(outputDirName + "/final-results/log.txt")));
                    Process process = builder.start();
                    BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
                    while (reader.readLine() != null) {
                    }
                    process.waitFor();
                } catch (Exception e) {
                    e.printStackTrace();
                }
                //for each genome in bin, get fasta
                for (int genomeId:genomesInBin) 
                {
                    cmd = "";
                    try 
                    {
                        cmd = "cd " + outputDirName + "/final-results\n"
                                + "cd scaffold-" + sortedScaffoldNum + "\n"
                                + "bash filterbyname.sh "
                                + "in=" + dbDir + "/" + dbFastaFile
                                + " out=" + genomeId + ".fasta names=" + indexToContigMap.get(genomeId)
                                + " include=t\n";
                        
                        FileWriter shellFileWriter = new FileWriter(outputDirName + "/final-results/run.sh");
                        shellFileWriter.write("#!/bin/bash\n");
                        shellFileWriter.write(cmd);
                        shellFileWriter.close();

                        ProcessBuilder builder = new ProcessBuilder("sh", outputDirName + "/final-results/run.sh");
                        builder.redirectError(Redirect.appendTo(new File(outputDirName + "/final-results/log.txt")));
                        Process process = builder.start();
                        BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
                        while (reader.readLine() != null) {
                        }
                        process.waitFor();
                    } catch (Exception ex) {
                        ex.printStackTrace();
                    }
                }   
                //for each genome in bin, get ANI
                for (int genomeId:genomesInBin) 
                {
                    ANIRes = "";
                    File genomeFastaFile = new File(outputDirName + "/final-results/scaffold-" + sortedScaffoldNum
                            + "/" + genomeId + ".fasta");
                    if(!genomeFastaFile.exists()) {
                        System.out.println("Could not find genome fasta file");
                    }
                    cmd = "";
                    
                    try {
                        cmd = "cd " + outputDirName + "/final-results\n"
                                + "cd scaffold-" + sortedScaffoldNum + "\n"
                                + "rm dnadiff-output.*\n"
                                + "dnadiff "
                                + genomeId 
                                + ".fasta scaffold-" + sortedScaffoldNum + ".fasta -p dnadiff-output\n";
                        
                        FileWriter shellFileWriter = new FileWriter(outputDirName + "/final-results/run.sh");
                        shellFileWriter.write("#!/bin/bash\n");
                        shellFileWriter.write(cmd);
                        shellFileWriter.close();
        
                        ProcessBuilder builder = new ProcessBuilder("sh", outputDirName + "/final-results/run.sh");
                        builder.redirectError(Redirect.appendTo(new File(outputDirName + "/final-results/log.txt")));
                        Process process = builder.start();
                        BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
                        while (reader.readLine() != null) {
                        }
                        process.waitFor();
                    } catch (Exception e) {
                        e.printStackTrace();
                    }
                    
                    refANI = 0;
                    qryANI = 0;
                    alignedNuclRef = 0;
                    alignedNuclQry = 0;
                    refLength = 0;
                    qryLength = 0;
                    
                    try {
                        brDnaDiff = new BufferedReader(new FileReader(outputDirName + "/final-results/scaffold-" 
                                + sortedScaffoldNum + "/dnadiff-output.report"));
                        while ((str = brDnaDiff.readLine()) != null) {
                            str = str.trim();
                            if (str.startsWith("TotalBases")) {
                                results = str.split("\\s+");
                                refLength = Integer.parseInt(results[1]);
                                qryLength = Integer.parseInt(results[2]);
                            }
                            else if (str.startsWith("AlignedBases")) {
                                results = str.split("\\s+");
                                String alignedRef = results[1].substring(results[1].indexOf('(')+1, results[1].length()-2);
                                alignedNuclRef = Double.parseDouble(alignedRef);
                                String alignedQry = results[2].substring(results[2].indexOf('(')+1, results[2].length()-2);
                                alignedNuclQry = Double.parseDouble(alignedQry);
                            }
                            else if (str.startsWith("AvgIdentity")) {
                                results = str.split("\\s+");
                                refANI = Double.parseDouble(results[1]);
                                qryANI = Double.parseDouble(results[2]);
                                break;
                            }
                        }
                        
                        ANIRes = "Reference:" + indexToContigMap.get(genomeId) + " Length: " + refLength
                                + " bp (ANI:" + refANI
                                + "%, Aligned Nucleotide:" + alignedNuclRef + "%), Query:scaffold-" 
                                + sortedScaffoldNum + ".fasta" + " Length: " + qryLength 
                                + " bp (ANI:" + qryANI
                                + "%, Aligned Nucleotide:" + alignedNuclQry + "%)";
                        bwANI.write("Scaffold-" + sortedScaffoldNum + "\t" + fastaHeader + "\t" + ANIRes + "\n");
                        brDnaDiff.close();
                    } catch (Exception e) {
                        e.printStackTrace();
                    }
                }               
                sortedScaffoldNum++;
            }
            
            br.close();
            bwBin.close();
            bwANI.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("Total final scaffolds: " + sortedScaffoldNum);
    }
	
	private void generateCoverageStat() {
		String cmd = "";
		
		Map<Integer, Integer> scaffoldToBinMap = new HashMap<Integer, Integer>();
		//get scaffold-bin map
		BufferedReader br = null;
		try {
			br = new BufferedReader(new FileReader(outputDirName + "/seed-scaffolds/binInfo.txt"));
			String str = "";
			String[] results;
			while ((str = br.readLine()) != null) 
			{
				str = str.trim();
				results = str.split("\t");
				int scaffold = Integer.parseInt(results[0].substring(results[0].indexOf('-')+1));
				int bin = Integer.parseInt(results[2].substring(results[2].indexOf('-')+1));
				scaffoldToBinMap.put(scaffold, bin);
			}
			br.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		//System.out.println("Scaffold to bin map size: " + scaffoldToBinMap.size());
		
		/*int sortedScaffoldNum = 0;
		try {
			br = new BufferedReader(new FileReader(outputDirName 
					+ "/final-results/final-scaffolds.fasta.fai"));
			String str = "";
			
			while ((str = br.readLine()) != null) 
			{
				str = str.trim();
				
				cmd = "";
				
				try {
					cmd = "cd " + outputDirName + "/final-results\n"
							+ "cd scaffold-" + sortedScaffoldNum + "\n"
							+ "awk '{print $1\"\tN/A\tN/A\t\"$2}' scaffold-" + sortedScaffoldNum + ".fasta.fai > length-list.txt\n";
					
					FileWriter shellFileWriter = new FileWriter(outputDirName + "/final-results/run.sh");
					shellFileWriter.write("#!/bin/bash\n");
					shellFileWriter.write(cmd);
					shellFileWriter.close();
	
					ProcessBuilder builder = new ProcessBuilder("sh", outputDirName + "/final-results/run.sh");
					builder.redirectError(Redirect.appendTo(new File(outputDirName + "/final-results/log.txt")));
					Process process = builder.start();
					BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
					while (reader.readLine() != null) {
					}
					process.waitFor();
				} catch (Exception e) {
					e.printStackTrace();
				}				
				sortedScaffoldNum++;
			}
			br.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		//System.out.println("Total scaffolds: " + sortedScaffoldNum);*/
		
		int sortedScaffoldNum = 0;
        try {
            br = new BufferedReader(new FileReader(outputDirName 
                    + "/final-results/final-scaffolds.fasta.fai"));
            String str = "";
            
            while ((str = br.readLine()) != null) 
            {
                sortedScaffoldNum++;
            }
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        
		ExecutorService executor = Executors.newFixedThreadPool(numThreads);
		for (int i = 0; i < sortedScaffoldNum; i++) {
			Runnable worker = new MyRunnableGetCoverageInfo(i);
			executor.execute(worker);
		}
		executor.shutdown();
		//wait until all threads are finished
		while (!executor.isTerminated()) {
			
		}
		//System.out.println("Finished all threads");
	}
	
	public class MyRunnableGetCoverageInfo implements Runnable {
		private final int sortedScaffoldNum;
		
		public MyRunnableGetCoverageInfo(int sortedScaffoldNum) {
			this.sortedScaffoldNum = sortedScaffoldNum;
		}
		
		@Override
		public void run() {
			getCoverageInfo(sortedScaffoldNum);
		}
	}
	
	private void getCoverageInfo(int sortedScaffoldNum) {
		String cmd = "";
		
		try {
		    if (read2.isEmpty()) {
		        cmd = "cd " + outputDirName + "/final-results\n"
                    + "cd scaffold-" + sortedScaffoldNum + "\n"
                    + "bowtie2-build scaffold-" + sortedScaffoldNum + ".fasta bowtie2-index\n"
                    + "bowtie2 -x bowtie2-index "
                    + "-U " + read1
                    + " | samtools view -bS - | samtools view -h -F 0x04 -b - | "
                    + "samtools sort - -o bowtie2-mapped.sam\n"
                    + "samtools view -bS bowtie2-mapped.sam | samtools sort - -o bowtie2-mapped.bam\n"
                    + "samtools depth bowtie2-mapped.bam > samtools-coverage.txt\n"
                    + "bash pileup.sh "
                    + "in=bowtie2-mapped.sam out=bbmap-coverage.txt\n";
		    }
		    else {
        			cmd = "cd " + outputDirName + "/final-results\n"
                    + "cd scaffold-" + sortedScaffoldNum + "\n"
                    + "bowtie2-build scaffold-" + sortedScaffoldNum + ".fasta bowtie2-index\n"
                    + "bowtie2 -x bowtie2-index "
                    + "-1 " + read1
                    + " -2 " + read2
                    + " | samtools view -bS - | samtools view -h -F 0x04 -b - | "
                    + "samtools sort - -o bowtie2-mapped.sam\n"
                    + "samtools view -bS bowtie2-mapped.sam | samtools sort - -o bowtie2-mapped.bam\n"
                    + "samtools depth bowtie2-mapped.bam > samtools-coverage.txt\n"
                    + "bash pileup.sh "
                    + "in=bowtie2-mapped.sam out=bbmap-coverage.txt\n";
		    }
			
			FileWriter shellFileWriter = new FileWriter(outputDirName +
					"/final-results/scaffold-" + sortedScaffoldNum + "/run.sh");
			shellFileWriter.write("#!/bin/bash\n");
			shellFileWriter.write(cmd);
			shellFileWriter.close();

			ProcessBuilder builder = new ProcessBuilder("sh", outputDirName +
					"/final-results/scaffold-" + sortedScaffoldNum + "/run.sh");
			builder.redirectError(Redirect.appendTo(new File(outputDirName +
					"/final-results/scaffold-" + sortedScaffoldNum + "/log.txt")));
			Process process = builder.start();
			BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
			while (reader.readLine() != null) {
			}
			process.waitFor();
		} catch (Exception e) {
			e.printStackTrace();
		}
		//System.out.println("Finished processing " + sortedScaffoldNum);
	}
	
	private void summerizeCoverageInfo() {
		BufferedReader br1 = null;
		BufferedReader br2 = null;
		
		BufferedWriter bw = null;
		try {
			bw = new BufferedWriter(new FileWriter(outputDirName 
					+ "/final-results/coverage-info.tsv"));
			
			br1 = new BufferedReader(new FileReader(outputDirName 
					+ "/final-results/bin-info.tsv"));
			
			br1.readLine();
			String str = "";
			String[] results;
			int scaffoldNum = 0;
			String scaffInfo = "";
			String bbmapCov = "";
			br1.readLine();
			
			while ((str = br1.readLine()) != null) 
			{
				str = str.trim();
				results = str.split("\t");
				scaffInfo = results[0].trim() + results[1].trim();
				try {
					br2 = new BufferedReader(new FileReader(outputDirName 
							+ "/final-results/scaffold-" + scaffoldNum + "/bbmap-coverage.txt"));
					br2.readLine();
					bbmapCov = br2.readLine();
					br2.close();
				} catch (Exception e) {
					e.printStackTrace();
				}
				
				bw.write(scaffInfo + "\t" + bbmapCov + "\n");
				scaffoldNum++;
			}
			br1.close();
			bw.close();
			System.out.println("Total final scaffolds: " + scaffoldNum);
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}	
	
	public void getANIInfo() {
		getBinAndANIInfo();
	}
	
	public void getCoverageInfo() {
		generateCoverageStat();
		summerizeCoverageInfo();
	}
	
	public void initialize(String read1, String read2, String outputDirName, 
			String FVEResDir, String dbType, String dbDir, int numThreads, int minFoldCov, int minScaffoldLen,
			int topBins, String spadesKmerLength) {
		
		this.read1 = read1;
		this.read2  = read2;
		this.outputDirName = outputDirName;
		this.FVEResDir = FVEResDir;
		this.dbType = dbType;
		this.dbDir = dbDir;
		this.numThreads = numThreads;
		this.minFoldCov = minFoldCov;
		this.minScaffoldLen = minScaffoldLen;
		this.topBins = topBins;
		this.spadesKmerlen = spadesKmerLength;
		
		setDbFiles();
		generateCoverageStatForFVE();
		getReadLen();
		createHashMap();
		readFVERes();
		getClusters();
		checkCoverage();
		getTopBins();
	}
}
