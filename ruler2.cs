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

///
public class netrule
{
  ///
  private class tf
  {
	  public float f1;
	  public float f2;
	  public float f3;
  }

  private class threeint
  {
  public int i1;
  public int i2;
  public int i3;
	public threeint()
	{
	}
	public threeint(int ai1, int ai2, int ai3)
	{
		i1 = ai1;
		i2 = ai2;
		i3 = ai3;
	}
  }


  ///
  private int quality;
  ///
  private string name;
  ///
  private Array<Point2d> points = new Array<Point2d>();
  ///
  private Array<INDEX_2> lines = new Array<INDEX_2>();
  ///
  private Array<Point2d> freezone = new Array<Point2d>();
  private Array<Point2d> freezonelimit = new Array<Point2d>();
  ///
  private Array<Array<Point2d>> freezone_i = new Array<Array<Point2d>>();
  ///
  private Array<Point2d> transfreezone = new Array<Point2d>();

  ///
  private Array<int> dellines = new Array<int>();
  ///
  private Array<Element2d> elements = new Array<Element2d>();
  ///
  private Array<tf> tolerances = new Array<tf>();
  private Array<tf> linetolerances = new Array<tf>();
  ///
  private Array<threeint> orientations = new Array<threeint>();
  ///
  private DenseMatrix oldutonewu = new DenseMatrix();
  private DenseMatrix oldutofreearea = new DenseMatrix();
  private DenseMatrix oldutofreearealimit = new DenseMatrix();
  ///
  private Array<DenseMatrix> oldutofreearea_i = new Array<DenseMatrix>();
  ///
  private MatrixFixWidth < 3> freesetinequ;

  ///
  private Array<Vec2d> linevecs = new Array<Vec2d>();

  ///
  private int noldp;
  private int noldl;
  ///
  private float fzminx;
  private float fzmaxx;
  private float fzminy;
  private float fzmaxy;

  /// topological distance of line to base element
  private Array<int> lnearness = new Array<int>();


  ///
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  netrule();
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
//ORIGINAL LINE: int GetNL() const
  public int GetNL()
  {
	  return lines.Size();
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
//ORIGINAL LINE: int GetNOldP() const
  public int GetNOldP()
  {
	  return noldp;
  }
  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetNOldL() const
  public int GetNOldL()
  {
	  return noldl;
  }
  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetNDelL() const
  public int GetNDelL()
  {
	  return dellines.Size();
  }
  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetNOrientations() const
  public int GetNOrientations()
  {
	  return orientations.Size();
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
//ORIGINAL LINE: int GetLNearness(int li) const
  public int GetLNearness(int li)
  {
	  return lnearness.Get(li);
  }

  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const Point2d & GetPoint(int i) const
  public Point2d GetPoint(int i)
  {
	  return points.Get(i);
  }
  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const INDEX_2 & GetLine(int i) const
  public INDEX_2 GetLine(int i)
  {
	  return lines.Get(i);
  }
  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const Element2d & GetElement(int i) const
  public Element2d GetElement(int i)
  {
	  return elements.Get(i);
  }
  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const threeint & GetOrientation(int i) const
  public threeint GetOrientation(int i)
  {
	  return orientations.Get(i);
  }
  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetDelLine(int i) const
  public int GetDelLine(int i)
  {
	  return dellines.Get(i);
  }
  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const Array<int> & GetDelLines() const
  public Array<int> GetDelLines()
  {
	  return dellines;
  }
  ///
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  void GetFreeZone(Array<Point2d> afreearea);
  ///

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: double CalcPointDist(int pi, const Point2d & p) const
  public double CalcPointDist(int pi, Point2d p)
  {
	double dx = p.X() - points.Get(pi).X();
	double dy = p.Y() - points.Get(pi).Y();
	tf tfp = tolerances.Get(pi);
	return tfp.f1 * dx * dx + tfp.f2 * dx * dy + tfp.f3 * dy * dy;
  }

  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: float CalcLineError(int li, const Vec2d & v) const;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  float CalcLineError(int li, Vec2d v);

  ///
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  void SetFreeZoneTransformation(Vector u, int tolclass);

  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool IsInFreeZone(const Point2d & p) const
  public bool IsInFreeZone(Point2d p)
  {
	if (p.X() < fzminx || p.X() > fzmaxx || p.Y() < fzminy || p.Y() > fzmaxy)
	{
		return false;
	}

	for (int i = 0; i < transfreezone.Size(); i++)
	{
	if (freesetinequ(i, 0) * p.X() + freesetinequ(i, 1) * p.Y() + freesetinequ(i, 2) > 0)
	{
		return false;
	}
	}
	return true;
  }

  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int IsLineInFreeZone(const Point2d & p1, const Point2d & p2) const
  public int IsLineInFreeZone(Point2d p1, Point2d p2)
  {
	if ((p1.X() > fzmaxx && p2.X() > fzmaxx) || (p1.X() < fzminx && p2.X() < fzminx) || (p1.Y() > fzmaxy && p2.Y() > fzmaxy) || (p1.Y() < fzminy && p2.Y() < fzminy))
	{
		return 0;
	}
	return IsLineInFreeZone2(p1, p2);
  }
  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int IsLineInFreeZone2(const Point2d & p1, const Point2d & p2) const;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  int IsLineInFreeZone2(Point2d p1, Point2d p2);
  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int ConvexFreeZone() const;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  int ConvexFreeZone();
  ///
  public Array<Point2d> GetTransFreeZone()
  {
	  return transfreezone;
  }

  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetPointNr(int ln, int endp) const
  public int GetPointNr(int ln, int endp)
  {
	  return lines.Get(ln).I(endp);
  }

  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const DenseMatrix & GetOldUToNewU() const
  public DenseMatrix GetOldUToNewU()
  {
	  return oldutonewu;
  }
  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const DenseMatrix & GetOldUToFreeArea() const
  public DenseMatrix GetOldUToFreeArea()
  {
	  return oldutofreearea;
  }
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
}



namespace netgen
{



//C++ TO C# CONVERTER NOTE: Enums must be named in C#, so the following enum has been named by the converter:
  public enum AnonymousEnum
  {
	  MAX_NEARNESS = 3
  }






}
