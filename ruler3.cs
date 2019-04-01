//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define DLL_HEADER __declspec(dllexport)
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define DLL_HEADER __declspec(dllimport)
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define __assume(cond) if (!(cond)) __builtin_unreachable(); else;
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define __assume(cond)
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define NG_INLINE __forceinline inline
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define NG_INLINE __forceinline inline
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define NG_INLINE inline
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define NG_INLINE inline
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define NGCORE_API_EXPORT __declspec(dllexport)
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define NGCORE_API_IMPORT __declspec(dllimport)
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define NGCORE_API NGCORE_API_EXPORT
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define NGCORE_API NGCORE_API_IMPORT
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define NETGEN_INLINE __forceinline inline
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define NETGEN_INLINE __forceinline inline
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define NETGEN_INLINE inline
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define NETGEN_INLINE inline
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define NETGEN_CORE_NGEXEPTION_STR_HELPER(x) #x
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define NETGEN_CORE_NGEXEPTION_STR(x) NETGEN_CORE_NGEXEPTION_STR_HELPER(x)
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define NG_EXCEPTION(s) ngcore::Exception(__FILE__ ":" NETGEN_CORE_NGEXEPTION_STR(__LINE__) "\t"+std::string(s))
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define NETGEN_CHECK_RANGE(value, min, max) { if ((value)<(min) || (value)>=(max)) throw ngcore::RangeException(__FILE__ ":" NETGEN_CORE_NGEXEPTION_STR(__LINE__) "\t", (value), (min), (max)); }
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define NETGEN_CHECK_RANGE(value, min, max)
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define NETGEN_DEBUG_LOG(logger, ...) SPDLOG_DEBUG(logger, __VA_ARGS__)
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define NETGEN_DEBUG_LOG(logger, ...)
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define VT_USER_START(n)
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define VT_USER_END(n)
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define VT_TRACER(n)


/**
  3D element generation rule.
 */
public class vnetrule
{
  /// rule is applicable for quality classes above this value
  private int quality;
  /// name of rule
  private string name;
  /// point coordinates in reference position
  private Array<Point3d> points = new Array<Point3d>();
  /// old and new faces in reference numbering
  private Array<Element2d> faces = new Array<Element2d>();
  /// additional edges of rule
  private Array<twoint> edges = new Array<twoint>();

  /// points of freezone in reference coordinates
  private Array<Point3d> freezone = new Array<Point3d>();
  /// points of freezone in reference coordinates if tolcalss to infty
  private Array<Point3d> freezonelimit = new Array<Point3d>();
  /// point index, if point equal to mappoint, otherwise 0
  private Array<int> freezonepi = new Array<int>();
  /// faces of each convex part of freezone
  private Array<Array<threeint>> freefaces = new Array<Array<threeint>>();
  /// set of points of each convex part of freezone
  private Array<Array<int>> freesets = new Array<Array<int>>();
  /// points of transformed freezone
  private Array<Point3d> transfreezone = new Array<Point3d>();
  /// edges of each convex part of freezone
  private Array<Array<twoint>> freeedges = new Array<Array<twoint>>();

  /// face numbers to be deleted
  private Array<int> delfaces = new Array<int>();
  /// elements to be generated
  private Array<Element> elements = new Array<Element>();
  /// tolerances for points and faces (used ??)
  private Array<double> tolerances = new Array<double>();
  private Array<double> linetolerances = new Array<double>();
  /// transformation matrix 
  private DenseMatrix oldutonewu = new DenseMatrix();
  /// transformation matrix: deviation old point to dev. freezone
  private DenseMatrix oldutofreezone;
  /** transformation matrix: deviation old point to dev. freezone, 
    quality class to infinity */
  private DenseMatrix oldutofreezonelimit;

  // can be deleted:
  // BaseMatrix *outf, *outfl;

  /**
    a point is outside of convex part of freezone, 
    iff mat * (point, 1) >= 0 for each component (correct ?)
    */
  private Array<DenseMatrix> freefaceinequ = new Array<DenseMatrix>();
  /// 
  private Array<fourint> orientations = new Array<fourint>();
  /**
    flags specified in rule-description file:
    t .. test rule
    */
  private Array<char> flags = new Array<char>();

  /**
    topological distance of face to base element
    non-connected: > 100  (??) 
    */
  private Array<int> fnearness = new Array<int>();
  private Array<int> pnearness = new Array<int>();
  private int maxpnearness;

  /// number of old points in rule
  private int noldp;
  /// number of new poitns in rule
  private int noldf;
  /// box containing free-zone
  // double fzminx, fzmaxx, fzminy, fzmaxy, fzminz, fzmaxz;
  public Box3d fzbox = new Box3d();


  ///
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  vnetrule();
  ///
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  public void Dispose();
  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetNP() const
  public int GetNP()
  {
	  return points.Size();
  }
  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetNF() const
  public int GetNF()
  {
	  return faces.Size();
  }
  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetNE() const
  public int GetNE()
  {
	  return elements.Size();
  }
  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetNO() const
  public int GetNO()
  {
	  return orientations.Size();
  }
  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetNEd() const
  public int GetNEd()
  {
	  return edges.Size();
  }
  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetNOldP() const
  public int GetNOldP()
  {
	  return noldp;
  }
  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetNOldF() const
  public int GetNOldF()
  {
	  return noldf;
  }
  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetNDelF() const
  public int GetNDelF()
  {
	  return delfaces.Size();
  }
  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetQuality() const
  public int GetQuality()
  {
	  return quality;
  }
  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetFNearness(int fi) const
  public int GetFNearness(int fi)
  {
	  return fnearness.Get(fi);
  }
  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetPNearness(int pi) const
  public int GetPNearness(int pi)
  {
	  return pnearness.Get(pi);
  }
  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetMaxPNearness() const
  public int GetMaxPNearness()
  {
	  return maxpnearness;
  }


  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const Point3d & GetPoint(int i) const
  public Point3d GetPoint(int i)
  {
	  return points.Get(i);
  }
  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const Element2d & GetFace(int i) const
  public Element2d GetFace(int i)
  {
	  return faces.Get(i);
  }
  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const Element & GetElement(int i) const
  public Element GetElement(int i)
  {
	  return elements.Get(i);
  }
  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const twoint & GetEdge(int i) const
  public twoint GetEdge(int i)
  {
	  return edges.Get(i);
  }
  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetDelFace(int i) const
  public int GetDelFace(int i)
  {
	  return delfaces.Get(i);
  }
  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int IsDelFace(int fn) const;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  int IsDelFace(int fn);

  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: float CalcPointDist(int pi, const Point3d & p) const;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  float CalcPointDist(int pi, Point3d p);
  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: double PointDistFactor(int pi) const
  public double PointDistFactor(int pi)
  {
	  return tolerances.Get(pi);
  }
  ///
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  void SetFreeZoneTransformation(Vector allp, int tolclass);
  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int IsInFreeZone(const Point3d & p) const;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  int IsInFreeZone(Point3d p);
  /**
    0 not in free-zone
    1 in free-zone
    -1 maybe 
   */
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  int IsTriangleInFreeZone(Point3d p1, Point3d p2, Point3d p3, Array<int> pi, int newone);
  ///
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  int IsQuadInFreeZone(Point3d p1, Point3d p2, Point3d p3, Point3d p4, Array<int> pi, int newone);
  ///
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  int IsTriangleInFreeSet(Point3d p1, Point3d p2, Point3d p3, int fs, Array<int> pi, int newone);

  ///
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  int IsQuadInFreeSet(Point3d p1, Point3d p2, Point3d p3, Point3d p4, int fs, Array<int> pi, int newone);

  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int ConvexFreeZone() const;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  int ConvexFreeZone();

  /// if t1 and t2 are neighbourtriangles, NTP returns the opposite Point of t1 in t2
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int NeighbourTrianglePoint(const threeint & t1, const threeint & t2) const;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  int NeighbourTrianglePoint(threeint t1, threeint t2);
  ///
  public Point3d GetTransFreeZone(int i)
  {
	  return transfreezone.Get(i);
  }

  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetNP(int fn) const
  public int GetNP(int fn)
  {
	  return faces.Get(fn).GetNP();
  }
  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetPointNr(int fn, int endp) const
  public int GetPointNr(int fn, int endp)
  {
	  return faces.Get(fn).PNum(endp);
  }
  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetPointNrMod(int fn, int endp) const
  public int GetPointNrMod(int fn, int endp)
  {
	  return faces.Get(fn).PNumMod(endp);
  }
  ///
  public fourint GetOrientation(int i)
  {
	  return orientations.Get(i);
  }

  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int TestFlag(char flag) const;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  int TestFlag(char flag);

  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const DenseMatrix & GetOldUToNewU() const
  public DenseMatrix GetOldUToNewU()
  {
	  return oldutonewu;
  }
  //
  //  const DenseMatrix & GetOldUToFreeZone () const { return oldutofreezone; }
  //
  //  const DenseMatrix & GetOldUToFreeZoneLimit () const 
  //    { return oldutofreezonelimit; }
  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const char * Name() const
  public string Name()
  {
	  return name;
  }
  ///
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  void LoadRule(istream ist);

  ///
  public Array<Point3d> GetTransFreeZone()
  {
	  return transfreezone;
  }
  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int TestOk() const;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  int TestOk();

  ///
//C++ TO C# CONVERTER TODO TASK: C# has no concept of a 'friend' function:
//ORIGINAL LINE: friend void TestRules();
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  void TestRules();
  ///
  //  friend void Plot3DRule (const ROT3D & r, char key);
}







namespace netgen
{






//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
//C++ TO C# CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
Timer ApplyRules_t("ruler3 - all");
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
//C++ TO C# CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
Timer ApplyRules_tstart("ruler3 - rule start");
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
//C++ TO C# CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
Timer ApplyRules_tloop("ruler3 - rule loop");

}
