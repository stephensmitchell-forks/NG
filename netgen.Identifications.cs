namespace netgen
{

	public class Identifications
	{
//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public Identifications(Mesh amesh)
		  {
			  this.mesh = amesh;
			  this.identifiedpoints = 100;
			  this.identifiedpoints_nr = 100;
			// identifiedpoints = new INDEX_2_HASHTABLE<int>(100);
			// identifiedpoints_nr = new INDEX_3_HASHTABLE<int>(100);
			maxidentnr = 0;
		  }
	}
}