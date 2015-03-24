/**
 * Created by tiantian on 2/24/15.
 */
public class Exon {
    public String chr;
    public String name;
    public String type;
    public Long startCo;
    public Long endCo;
    public String score;
    public String strand;
    public String frameRead;
    public String attribute;

    public String toString() {
        return chr + "\t" + name + "\t" + type + "\t" + startCo + "\t" + endCo +
                "\t" + score + "\t" + strand + "\t" + frameRead + "\t" + attribute;
    }
}
