/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mergae.tauchiivariants.v02;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;

/**
 *
 * @author shahbazlums
 */
public class MergAeTauchiiVariantsV02 {

    /**
     * @param args the command line arguments
     * @throws java.io.IOException
     */
    public static void main(String[] args) throws IOException {
        // TODO code application logic here

        args = "merge -snp1 testdata/AE1055.samtools.vcf -bc1 testdata/AE1055.bc -snp2 testdata/AE1068.samtools.vcf -bc2 testdata/AE1068.bc -snp3 testdata/AE1650.samtools.vcf -bc3 testdata/AE1650.bc -out testdata/megeOutput.txt".split(" +");
//        args = "merge -snp1 testdata/AE1055.samtools.vcf -snp2 testdata/AE1068.samtools.vcf -snp3 testdata/AE1650.samtools.vcf -bc1 testdata/AE1055.bc -bc2 testdata/AE1068.bc -bc3 testdata/AE1650.bc -out testdata/megeOutput.txt".split(" +");

        String SNPFiles = "";
        String CoverageFiles = "";
        String OutputFileName = "";
        int idx = 0;
        while(idx < args.length){
            
            if(args[idx].equals("merge")){
                idx++;
                //System.out.println(args[idx]);
            }else if(args[idx].startsWith("-snp")){
                SNPFiles = SNPFiles.concat(args[++idx] + " ");
            }else if(args[idx].startsWith("-bc")){
                CoverageFiles = CoverageFiles.concat(args[++idx] + " ");
            }else if(args[idx].startsWith("-out")){
                OutputFileName = args[++idx];
            }else{
                idx++;
            }
        }
        
// Part 1: Merge variants
        String[] SNPFileNames = SNPFiles.split(" ");
        String[] CoverageFileNames = CoverageFiles.split(" ");
        int CoverageThreshold = 3; 
        //int coverageTotal = 0;
        
        //Input Files check
        if(CoverageFiles.isEmpty() || SNPFiles.isEmpty()){
            printHelpMerge();
            return;
        }
        //Check
        if(CoverageFileNames.length != SNPFileNames.length){
            System.err.println("Either SNP or Coverage file is missing...");
            return;
        }
        //Check
        if(OutputFileName.isEmpty()){
            System.err.println("Please provide output filename!");
            return;
        }
        
//Write Header and SNPFilenames
        BufferedWriter writer = new BufferedWriter(new FileWriter(OutputFileName));
        writer.write("#Header\n#Scaffolds_id\tSNPpos\t");
        //System.out.print("#Header\n#Scaffolds_id\tSNPpos\t");
        
        
        SortedMap<String,Map<Integer, String[]>> scaffolds = new TreeMap();
        String line;        
//Read each VCF file and prepare a single merged VCF TreeMap
        for(int i=0;i<SNPFileNames.length;i++){
            String[] filename = SNPFileNames[i].split("/");                     //split the complete file path and just pick the filename
            String[] Accession = filename[filename.length-1].split("\\.|-");    //form that filename only pick accession
            writer.write(Accession[0] + "\t");                                  //write accession to output file
            //System.out.print(Accession[0]+"\t");        
        
            // read SNP file
            BufferedReader br = new BufferedReader(new FileReader(SNPFileNames[i]));
            while ((line = br.readLine()) != null) {                
                if(line.startsWith("#")){
                    continue;
                }

                String entries[] = line.split("\t");
                if(Float.parseFloat(entries[5]) < 20){ //if any SNP phred score is less than 20 then skip it
                    continue;
                }else if(entries[8].startsWith("GT") && entries[9].startsWith("0/1")){
                    entries[4] = "*";
                }
                //else if(entries[7].startsWith("INDEL")){
                //    continue;
                //}

                String[] array = new String[SNPFileNames.length];
                array[i] = entries[4];    
    
                //add variants with respect to their scaffolds
                if(!scaffolds.containsKey(entries[0])){
                    Map<Integer, String[]> variant1 = new TreeMap();
                    variant1.put(Integer.parseInt(entries[1]), array);
                    scaffolds.put(entries[0], variant1);
                }else{                    
                    Map<Integer, String[]> variants = scaffolds.get(entries[0]);
                    Integer varPos = Integer.parseInt(entries[1]);
                    if(variants.containsKey(varPos)){                        
                        String[] modifiedArr = variants.get(varPos);
                        modifiedArr[i] = entries[4];
                        variants.put(Integer.parseInt(entries[1]), modifiedArr);
                    }else{
                        variants.put(Integer.parseInt(entries[1]), array);                        
                    }                    
                    scaffolds.put(entries[0], variants);             
                    
                }                
            }//end of while loop
             
        }//end of for loop
        //System.out.println("");
        writer.write("\n");
        
//Part 2:process Coverage files and assign nucleotides to missing SNPs
        for(int i=0;i<CoverageFileNames.length;i++){//Pick coverage files one by one
        String Line;        
        BufferedReader coverageFileReader = new BufferedReader(new FileReader(CoverageFileNames[i]));
            while ((Line = coverageFileReader.readLine()) != null) {//Process coverage file
                if(Line.startsWith("@")){
                    continue;
                }  
                String[] tmpArray = Line.split("\t");
                String seqName = tmpArray[0];
                int position = Integer.parseInt(tmpArray[1]);
                
                //skip scaffold which dont have any variant
                if(!scaffolds.containsKey(seqName)){
                    continue;
                }
                
                Map<Integer, String[]> variants = scaffolds.get(seqName);                                
                if(variants.containsKey(position)){
                    String[] modifiedArr2 = variants.get(position);
                    //now go for covergae
                    if(modifiedArr2[i] == null){    //if yes                      
                        Map<Character, Integer> SNP_Coverage = new TreeMap<>();
                        //process coverage for each variant and the get nucleotide whose coverage is maximum
                        for(int j = 2; j< tmpArray.length; j++){
                            //coverageTotal = coverageTotal + Integer.parseInt(tmpArray[j]);
                            switch (j) {                            
                                case 2:
                                    SNP_Coverage.put('A', Integer.parseInt(tmpArray[j]));
                                    break;
                                case 3:
                                    SNP_Coverage.put('T', Integer.parseInt(tmpArray[j]));
                                    break;
                                case 4:
                                    SNP_Coverage.put('G', Integer.parseInt(tmpArray[j]));
                                    break;
                                case 5:
                                    SNP_Coverage.put('C', Integer.parseInt(tmpArray[j]));
                                    break;
                                case 6:
                                    SNP_Coverage.put('N', Integer.parseInt(tmpArray[j]));
                                    break;
                                default:
                                    break;
                            }//end of switch case
                        }

                        //Pick max covarage                    
                        Map.Entry<Character, Integer> maxEntry = null;
                        for (Map.Entry<Character, Integer> entry2 : SNP_Coverage.entrySet()){
                            if (maxEntry == null || entry2.getValue().compareTo(maxEntry.getValue()) > 0){
                                maxEntry = entry2;
                            }
                        }
                        
                        //Check for threshold and update the MergedSNP list
                        if(maxEntry.getValue() >= CoverageThreshold){
                            //System.out.println(maxEntry);                        
                            modifiedArr2[i] = Character.toString(maxEntry.getKey());
                            variants.put(position, modifiedArr2);          //update MergedSNPS
                        }else{
                            modifiedArr2[i] = "<";
                            variants.put(position, modifiedArr2);
                        }

                        //coverageTotal = 0;
                        
                    }
                }
                
                //System.out.println(seqName+"\t"+position);
              
            }//end of while loop
                        
        }//end of for loop 
        
        
        //System.out.println("");
//Part 3: Iterate over hashMap (merged SNPs)
        //First go for scaffold
        for (Map.Entry<String,Map<Integer, String[]>> scaffName : scaffolds.entrySet()) {
            //second fo for variant
            Map<Integer, String[]> variants_disp = scaffName.getValue();
            //thrid process each variant
            for (Integer key : variants_disp.keySet()){
                String[] arr = variants_disp.get(key); 
                //System.out.print(scaffName.getKey()+ "\t" + key);
                writer.write(scaffName.getKey()+ "\t" + key);
                //process all variants for single variant position
                for(String value : arr){
                    //System.out.println(entry.getKey()+ "\t"+ key + "\t" + value);
                    if(value == null){
                        value = "<";
                    }
                    //System.out.print("\t"+value);
                    writer.write("\t"+value);
                }
                //System.out.println("");
                writer.write("\n");
            }
            //System.out.println("");
            //System.out.println(entry.getKey()+ "\t"+out.keySet()+"\t"+ Arrays.toString(out.get(5494)));
        }
        writer.close();
    
    }
    private static void printHelpMerge() {
         System.out.println("Usage: java -jar MergeAllAe.tauchiiVariants.jar <command>\n");
         System.out.println("Command");
         System.out.println("\tmerge : Combine the Variants(vcf format) and uses coverage data to assign nucleotides to missing SNPs");
         System.out.println("\tExample: merge -snp1 diploid1.snp -snp2 diploid2.snp -snp3 diploid3.snp -bc1 diploid1.bc -bc2 diploid2.bc -bc3 diploid3.bc -out Output.vcf");
    }
}
