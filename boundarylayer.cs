/// Create a typical prismatic boundary layer on the given 
/// surfaces

public class BoundaryLayerParameters
{
  // parameters by Philippose ..
  public Array<int> surfid = new Array<int>();
  public Array<double> heights = new Array<double>();
  public Array<double> new_matnrs = new Array<double>();
  public int prismlayers = 1;
  public int bulk_matnr = 1;
  public int new_matnr = 1;
  public double hfirst = 0.01;
  public double growthfactor = 1;
  public bool optimize = true;
}





