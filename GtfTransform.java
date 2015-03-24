/**
 * Created by tiantian on 2/24/15.
 * This is a program that reads data from refSeq gene table and transform the data to a gtf file.
 * In this file, all isoforms of one gene is combined by union their exons.
 */
import java.util.*;
import java.io.*;

public class GtfTransform {
    public static void main(String[] args) throws FileNotFoundException {
        if(args.length != 2) {
            System.out.println("Usage: java GtfTransform inputName outputName");
            System.exit(0);
        }
        File geneData = new File(args[0]);
        Scanner input = new Scanner(geneData);

        input.nextLine();       // Discard the header line

        HashMap<Integer, Gene> genome = new HashMap<Integer, Gene>();   // Here key is the hash value of the string of gene symbol

        while(true) {
            String line = input.nextLine();     // read next line
            String[] row = line.split("\\t");   // split line by tab

            // If prefix of transcriptID is "NR" or "XR", means this is noncoding gene, then skip
            String transcriptID = row[1];
            if(transcriptID.substring(0,2).equals("NR") || transcriptID.substring(0,2).equals("XR")) {
                continue;
            }

            int exonNumber = Integer.parseInt(row[8]);      // get exon number
            String[] startCoos = row[9].split(",");         // get the array of exon start coordinates
            String[] endCoos = row[10].split(",");          // get the array of exon end coordinates

            if(exonNumber != startCoos.length || exonNumber != endCoos.length) {
                System.out.println(line);
                throw new AssertionError();
            }

            for(int i = 0; i<exonNumber; i++) {
                Exon comingExon = new Exon();
                comingExon.chr = row[2];                                // get chromosome name
                comingExon.name = "hg19_refGene";                       // source feature from GTF field
                comingExon.type = "exon";                               // the name of type of the feature
                comingExon.startCo = Long.parseLong(startCoos[i]);      // start coordinate
                comingExon.endCo = Long.parseLong(endCoos[i]);          // end coordinate
                comingExon.strand = row[3];                             // strand
                comingExon.score = "0.0";                               // score
                // frame for coding exon. We don't need the frame information here, so set all fram to "."
                comingExon.frameRead = ".";
                comingExon.attribute = "gene_id \"" + row[12] + "\"; transcript_id \"" + row[12] + "\";";

                // assert start coordinate is smaller than end coordinate
                if (comingExon.startCo.compareTo(comingExon.endCo) >= 0) {
                    System.out.println(comingExon);
                    throw new AssertionError();
                }

                // If the current exon belongs to a new gene
                // then create the gene and add it to the genome
                int symbolHash = row[12].hashCode();
                if (!genome.containsKey(symbolHash)) {
                    Gene newGene = new Gene();
                    newGene.exons.add(comingExon);
                    genome.put(symbolHash, newGene);
                }
                // Else add to the exist gene
                else {
                    Gene gene = genome.get(symbolHash);
                    gene.exons.add(comingExon);
                }
            }

            if(!input.hasNext()) break;
        }

        input.close();

        // merge exon by its coordinates
        for(Gene gene : genome.values()) {
            mergeExons(gene);
        }

        // sort gene by chromosome name
        ArrayList<Gene> genomeList = new ArrayList(genome.values());
        Collections.sort(genomeList, new Comparator<Gene>() {
           public int compare(Gene g1, Gene g2) {
               return g1.exons.get(0).chr.compareTo(g2.exons.get(0).chr);
           }
        });

        File outGtf = new File(args[1]);
        PrintWriter output = new PrintWriter(outGtf);

        for(Gene gene : genomeList) {
            output.write(gene.toString());
            output.flush();
        }
        output.close();
    }

    /**
     * Merges overlapping exons
     * Algorithm (reference: http://www.geeksforgeeks.org/merging-intervals/)
     * Here "intervals" means the region from start coordinate of exon to end coordinate
     * 1. Sort the intervals based on increasing order of starting time.
     * 2. Push the first interval on to a stack.
     * 3. For each interval do the following
     * ……..a. If the current interval does not overlap with the stack top, push it.
     * ……..b. If the current interval overlaps with stack top and ending time of current interval
     *        is more than that of stack top, update stack top with the ending time of current interval.
     * 4. At the end stack contains the merged intervals.
     * @param g
     */
    public static void mergeExons(Gene g) {
        // Test if the given gene has at least exon
        if(g.exons.size() <= 0) return;

        // Create an empty stack of intervals
        Stack<Exon> s = new Stack<Exon>();

        // sort the intervals based on start coordinate
        Collections.sort(g.exons, new Comparator<Exon>() {
            public int compare(Exon e1, Exon e2) {
                return e1.startCo.compareTo(e2.startCo);
            }
        });

        // push the first exon to stack
        s.push(g.exons.get(0));

        // Start from the next exon and merge if necessary
        for(int i=1; i<g.exons.size(); i++) {
            // get exon from stack top
            Exon top = s.peek();

            // if current exon is not overlapping with stack top,
            // push it to the stacl
            if(top.endCo < g.exons.get(i).startCo) {
                s.push(g.exons.get(i));
            }
            // Otherwise update the ending coordinate of top if ending of current
            // exon is more
            else if(top.endCo < g.exons.get(i).endCo) {
                top.endCo = g.exons.get(i).endCo;
                s.pop();
                s.push(top);
            }
        }

        // Reverse the order in stack
        Stack<Exon> reS = new Stack<Exon>();
        while(!s.empty()) {
            reS.push(s.pop());
        }

        // Put exons in stack to gene
        ArrayList<Exon> sortExon = new ArrayList<Exon>();
        while(!reS.empty()) {
            sortExon.add(reS.pop());
        }

        g.exons = sortExon;
    }
}
