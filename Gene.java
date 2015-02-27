/**
 * Created by tiantian on 2/24/15.
 */
import java.util.ArrayList;

public class Gene {
    public ArrayList<Exon> exons = new ArrayList<Exon>();

    public String toString() {
        String str = "";
        for(Exon e : exons) {
            str += e.toString() + "\n";
        }
        return str;
    }
}
