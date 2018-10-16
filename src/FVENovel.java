import java.time.LocalDateTime;

public class FVENovel {
	private static int MIN_FOLD_COV = 3;
	private static int MIN_SCAFFOLD_LEN = 2000;
	private static int TOP_BINS = 100;
	private static int MYTHREADS = 4;
	private static String read1 = "";
	private static String read2 = "";	
	private static String outputDirName = "";
	private static String FVEResDir = "";
	private static String dbType = "";
	private static String dbDir = "";
	private static boolean reportCoverage = false;
	private static String spadesKmerlen = "default";
	private static boolean generateSeeds = true;
	private static int seedExtensionStart = 0;
	private static boolean onlyStatistics = false;
	
	private static void printUsage() {
	    System.out.println("Usage:");
        System.out.println(
                "java -cp bin FVENovel -1 $read1File -2 $read2File -fveres $FastViromeExplorerResDirectory "
                + "-dbType gov -dbDir /path-to-referencedb-folder -o $outputDirectory");
        System.out.println("-1: input .fastq file for read sequences (paired-end 1), mandatory field.");
        System.out.println("-2: input .fastq file for read sequences (paired-end 2).");
        System.out.println("-fveres: full path of directory containing results from FastviromeExplorer.");
        System.out.println("-dbType: specify the database name (\"gov\", \"gov_epi_mes\", or \"imgvr\")");
        System.out.println("-dbDir: full path of directory containing all database files.");
        System.out.println("-o: full path of output directory.");
	}
	
	private static void parseArguments(String[] args) {
		if (args.length == 0) {
			printUsage();
			System.exit(1);
		} else {
			for (int i = 0; i < args.length; i++) {
				if (args[i].startsWith("-")) {
					if ((i + 1) >= args.length) {
						System.out.println("Missing argument after " + args[i] + ".");
						printUsage();
						System.exit(1);
					} else {
						if (args[i].equals("-o")) {
							outputDirName = args[i + 1];
						} else if (args[i].equals("-1")) {
							read1 = args[i + 1];
						} else if (args[i].equals("-2")) {
							read2 = args[i + 1];
						} else if (args[i].equals("-fveres")) {
							FVEResDir = args[i + 1];
						} else if (args[i].equals("-dbType")) {
							dbType = args[i + 1];
						} else if (args[i].equals("-dbDir")) {
							dbDir = args[i + 1];
						} else if (args[i].equals("-p")) {
							MYTHREADS = Integer.parseInt(args[i + 1]);
						} else if (args[i].equals("-minFoldCov")) {
							MIN_FOLD_COV = Integer.parseInt(args[i + 1]);
						} else if (args[i].equals("-minScaffoldLen")) {
							MIN_SCAFFOLD_LEN = Integer.parseInt(args[i + 1]);
						} else if (args[i].equals("-topBins")) {
							TOP_BINS = Integer.parseInt(args[i + 1]);
						} else if (args[i].equals("-spadeskmer")) {
                            spadesKmerlen = args[i + 1];
                        } else if (args[i].equals("-seedExtensionStart")) {
                            seedExtensionStart = Integer.parseInt(args[i + 1]);
                        } else if (args[i].equals("-generateSeeds")) {
                            if (args[i + 1].equals("true")) {
                                generateSeeds = true;
                            }
                            else if (args[i + 1].equals("false")) {
                                generateSeeds = false;
                            }
                            else {
                                generateSeeds = true;
                            }
                        } else if (args[i].equals("-reportCov")) {
							if (args[i + 1].equals("true")) {
								reportCoverage = true;
							}
							else {
								reportCoverage = false;
							}
						} else if (args[i].equals("-onlyStatistics")) {
                            if (args[i + 1].equals("true")) {
                                onlyStatistics = true;
                            }
                            else {
                                onlyStatistics = false;
                            }
                        } else {
							System.out.println("Invalid argument.");
							printUsage();
							System.exit(1);
						}
					}
				}
			}
		} // finish parsing arguments
	}	

	public static void main(String[] args) {
		parseArguments(args);
		System.out.println("Finished parsing command line arguments");
		
		GrowScaffolds growScaffolds = new GrowScaffolds();
		growScaffolds.initialize(read1, read2, outputDirName, FVEResDir, dbType, dbDir, MYTHREADS, MIN_FOLD_COV,
				MIN_SCAFFOLD_LEN, TOP_BINS, spadesKmerlen);
		if (onlyStatistics) {
		    growScaffolds.getANIInfo();
		    growScaffolds.getCoverageInfo();
		}
		else {
        		if (generateSeeds) {
                		System.out.println("Started 1st round of assembly: " + LocalDateTime.now());
                		growScaffolds.writeFastqFilesAndAssembly();
                		System.out.println("Finished 1st round of assembly: " + LocalDateTime.now());
        		}
        		else {
        		    System.out.println("Skipping 1st round of assembly");
        		}
        		
        		int totalSeedScaffolds = growScaffolds.getSeedScaffolds();
        		System.out.println("Total seed scaffolds: " + totalSeedScaffolds);
        		if (totalSeedScaffolds == 0) {
        			System.out.println("Total seed scaffolds: 0. Couldn't grow scaffolds.");
        			System.exit(0);
        		}
        		
        		System.out.println("Started growing seed scaffolds: " + LocalDateTime.now());
        		growScaffolds.growSeedScaffolds(totalSeedScaffolds, seedExtensionStart);
        		System.out.println("Finished growing seed scaffolds: " + LocalDateTime.now());
        		
        		growScaffolds.getFinalScaffolds(totalSeedScaffolds);
        		
        		growScaffolds.getANIInfo();
        		if (reportCoverage) {
        		    System.out.println("Started generating coverage stat for all scaffolds: " + LocalDateTime.now());
        			growScaffolds.getCoverageInfo();
        			System.out.println("Finished generating coverage stat for all scaffolds: " + LocalDateTime.now());
        		}
		}
	}
}
