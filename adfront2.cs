using System;

/*
  Advancing front class for surfaces
*/

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

/**************************************************************************/
/* File:   adfront2.hpp                                                   */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Okt. 95                                                    */
/**************************************************************************/


/**

    Advancing front class for surfaces

*/

  ///
  public class FrontPoint2
  {
	/// coordinates
	private Point < 3> p;
	/// global node index
	private PointIndex globalindex = new PointIndex();
	/// number of front lines connected to point 
	private int nlinetopoint;
	/// distance to original boundary
	private int frontnr;

	private bool onsurface;

	///
	public MultiPointGeomInfo[] mgi;

	///
	public FrontPoint2()
	{
	  globalindex.Invalidate(); //  = -1;
	  nlinetopoint = 0;
	  frontnr = INT_MAX - 10; // attention: overflow on calculating  INT_MAX + 1
	  mgi = null;
	  onsurface = true;
	}

	///
	public FrontPoint2(Point < 3> ap, PointIndex agi, MultiPointGeomInfo amgi, bool aonsurface = true)
	{
	  p = ap;
	  globalindex = agi;
	  nlinetopoint = 0;
	  frontnr = INT_MAX - 10;
	  onsurface = aonsurface;

	  if (amgi != null)
	  {
	  mgi = new MultiPointGeomInfo(amgi);
	  for (int i = 1; i <= mgi.GetNPGI(); i++)
	  {
		if (mgi.GetPGI(i).trignum <= 0)
		{
		  Console.Write("WARNING: Add FrontPoint2, illegal geominfo = ");
		  Console.Write(mgi.GetPGI(i).trignum);
		  Console.Write("\n");
		}
	  }
	  }
	  else
	  {
		mgi = null;
	  }
	}

	///
	public void Dispose()
	{
		;
	}

	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const Point<3> & P() const
	public Point < 3> P()
	{
		return p;
	}
	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: operator const Point<3> & () const
	public static implicit operator Point < 3> & (FrontPoint2 ImpliedObject)
	{
		return p;
	}
	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: PointIndex GlobalIndex() const
	public PointIndex GlobalIndex()
	{
		return globalindex;
	}

	///
	public void AddLine()
	{
		nlinetopoint++;
	}
	///
	public void RemoveLine()
	{
	  nlinetopoint--;
	  if (nlinetopoint == 0)
	  {
	nlinetopoint = -1;
	  }
	}

	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool Valid() const
	public bool Valid()
	{
		return nlinetopoint >= 0;
	}

	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool OnSurface() const
	public bool OnSurface()
	{
		return onsurface;
	}

	///
	public void DecFrontNr(int afrontnr)
	{
	  if (frontnr > afrontnr)
	  {
		  frontnr = afrontnr;
	  }
	}

	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int FrontNr() const
	public int FrontNr()
	{
		return frontnr;
	}
  }


  ///
  public class FrontLine
  {
	/// Point Indizes
	private INDEX_2 l = new INDEX_2();
	/// quality class 
	private int lineclass;
	/// geometry specific data
	private PointGeomInfo[] geominfo = Arrays.InitializeWithDefaultInstances<PointGeomInfo>(2);

	public FrontLine()
	{
	  lineclass = 1;
	}

	///
	public FrontLine(INDEX_2 al)
	{
	  l = al;
	  lineclass = 1;
	}


	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const INDEX_2 & L() const
	public INDEX_2 L()
	{
	  return l;
	}

	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int LineClass() const
	public int LineClass()
	{
	  return lineclass;
	}

	///
	public void IncrementClass()
	{
	  lineclass++;
	}
	///
	public void ResetClass()
	{
	  lineclass = 1;
	}

	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool Valid() const
	public bool Valid()
	{
	  return l.I1() != -1;
	}
	///
	public void Invalidate()
	{
	  l.I1() = -1;
	  l.I2() = -1;
	  lineclass = 1000;
	}

	public void SetGeomInfo(PointGeomInfo gi1, PointGeomInfo gi2)
	{
	geominfo[0] = gi1;
	geominfo[1] = gi2;
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const PointGeomInfo * GetGeomInfo() const
	public PointGeomInfo GetGeomInfo()
	{
		return geominfo;
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const PointGeomInfo & GetGeomInfo(int endp) const
	public PointGeomInfo GetGeomInfo(int endp)
	{
		return geominfo[endp - 1];
	}

//C++ TO C# CONVERTER TODO TASK: C# has no concept of a 'friend' class:
//	friend class AdFront2;
  }


public class AdFront2
{

  ///
  private Array<FrontPoint2> points = new Array<FrontPoint2>(); /// front points
  private Array<FrontLine> lines = new Array<FrontLine>(); /// front lines

  private Box3d boundingbox = new Box3d();
  private BoxTree < 3> linesearchtree; /// search tree for lines
  private Point3dTree pointsearchtree = new Point3dTree(); /// search tree for points
  private Point3dTree cpointsearchtree = new Point3dTree(); /// search tree for cone points (not used ???)

  private Array<int> delpointl = new Array<int>(); /// list of deleted front points
  private Array<int> dellinel = new Array<int>(); /// list of deleted front lines

  private int nfl; /// number of front lines;
  private INDEX_2_HASHTABLE<int> allflines; /// all front lines ever have been

  private Array<int> invpindex = new Array<int>();

  private int minval;
  private int starti;

  ///
  //  AdFront2 ();
  public AdFront2(Box3d aboundingbox)
  {
	  this.boundingbox = aboundingbox;
	  this.linesearchtree = new <type missing>(boundingbox.PMin(), boundingbox.PMax());
	  this.pointsearchtree = new <type missing>(boundingbox.PMin(), boundingbox.PMax());
	  this.cpointsearchtree = new <type missing>(boundingbox.PMin(), boundingbox.PMax());
	nfl = 0;
	allflines = 0;

	minval = 0;
	starti = lines.Begin();
  }

  ///
  public void Dispose()
  {
	if (allflines != null)
	{
		allflines.Dispose();
	}
  }

  ///
  // void GetPoints (Array<Point<3> > & apoints) const;
  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void Print(ostream & ost) const
  public void Print(ostream ost)
  {
	ost << points.Size() << " Points: " << "\n";
	for (int i = points.Begin(); i < points.End(); i++)
	{
	  if (points[i].Valid())
	  {
	ost << i << "  " << points[i].P() << "\n";
	  }
	}

	ost << nfl << " Lines: " << "\n";
	for (int i = lines.Begin(); i < lines.End(); i++)
	{
	  if (lines[i].Valid())
	  {
	ost << lines[i].L().I1() << " - " << lines[i].L().I2() << "\n";
	  }
	}

	ost << flush;
  }

  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool Empty() const
  public bool Empty()
  {
	return nfl == 0;
  }
  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetNFL() const
  public int GetNFL()
  {
	  return nfl;
  }

  public FrontLine GetLine(int nr)
  {
	  return lines[nr];
  }
  public FrontPoint2 GetPoint(int nr)
  {
	  return points[nr];
  }


  ///
  public int SelectBaseLine(ref Point < 3> p1, ref Point < 3> p2, PointGeomInfo geominfo1, PointGeomInfo geominfo2, ref int qualclass)
  {
	int baselineindex = -1;

	for (int i = starti; i < lines.End(); i++)
	{
	if (lines[i].Valid())
	{
		int hi = lines[i].LineClass() + points[lines[i].L().I1()].FrontNr() + points[lines[i].L().I2()].FrontNr();

		if (hi <= minval)
		{
		minval = hi;
		baselineindex = i;
		break;
		}
	}
	}

	if (baselineindex == -1)
	{
	minval = INT_MAX;
	for (int i = lines.Begin(); i < lines.End(); i++)
	{
	  if (lines[i].Valid())
	  {
		  int hi = lines[i].LineClass() + points[lines[i].L().I1()].FrontNr() + points[lines[i].L().I2()].FrontNr();

		  if (hi < minval)
		  {
		  minval = hi;
		  baselineindex = i;
		  }
	  }
	}
	}
	starti = baselineindex + 1;

	p1 = points[lines[baselineindex].L().I1()].P();
	p2 = points[lines[baselineindex].L().I2()].P();
	geominfo1 = lines[baselineindex].GetGeomInfo(1);
	geominfo2 = lines[baselineindex].GetGeomInfo(2);

	qualclass = lines[baselineindex].LineClass();

	return baselineindex;
  }

  ///
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private int GetLocals_timer = NgProfiler.CreateTimer("adfront2::GetLocals");

  public int GetLocals(int baselineindex, Array<Point3d> locpoints, Array<MultiPointGeomInfo> pgeominfo, Array<INDEX_2> loclines, Array<int> pindex, Array<int> lindex, double xh)
  {
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static int timer = NgProfiler::CreateTimer("adfront2::GetLocals");
	NgProfiler.RegionTimer reg = new NgProfiler.RegionTimer(GetLocals_timer);


	int pstind;
	Point < 3> midp, p0;

	pstind = lines[baselineindex].L().I1();
	p0 = points[pstind].P();

	loclines.Append(lines[baselineindex].L());
	lindex.Append(baselineindex);

	ArrayMem<int, 1000> nearlines = new ArrayMem<int, 1000>(0);
	ArrayMem<int, 1000> nearpoints = new ArrayMem<int, 1000>(0);

	// dominating costs !!
	linesearchtree.GetIntersecting(p0 - new Vec3d(xh, xh, xh), p0 + new Vec3d(xh, xh, xh), nearlines);

	pointsearchtree.GetIntersecting(p0 - new Vec3d(xh, xh, xh), p0 + new Vec3d(xh, xh, xh), nearpoints);

	for (int ii = 0; ii < nearlines.Size(); ii++)
	{
	int i = nearlines[ii];
	if (lines[i].Valid() && i != baselineindex)
	{
			loclines.Append(lines[i].L());
			lindex.Append(i);
	}
	}

	// static Array<int> invpindex;
	invpindex.SetSize(points.Size());
	// invpindex = -1;
	for (int i = 0; i < nearpoints.Size(); i++)
	{
	  invpindex[nearpoints[i]] = -1;
	}

	for (int i = 0; i < loclines.Size(); i++)
	{
	invpindex[loclines[i].I1()] = 0;
	invpindex[loclines[i].I2()] = 0;
	}


	for (int i = 0; i < loclines.Size(); i++)
	{
	for (int j = 0; j < 2; j++)
	{
		int pi = loclines[i][j];
		if (invpindex[pi] == 0)
		{
		pindex.Append(pi);
		invpindex[pi] = pindex.Size();
				locpoints.Append(points[pi].P());
		loclines[i][j] = locpoints.Size();
		}
		else
		{
		  loclines[i][j] = invpindex[pi];
		}
	}
	}


	// double xh2 = xh*xh;
	for (int ii = 0; ii < nearpoints.Size(); ii++)
	{
		int i = nearpoints[ii];
	if (points[i].Valid() && points[i].OnSurface() && invpindex[i] <= 0)
	{
			locpoints.Append(points[i].P());
		invpindex[i] = locpoints.Size();
		pindex.Append(i);
	}
	}
	/*
	double xh2 = xh*xh;
	for (i = 1; i <= points.Size(); i++)
	  {
	if (points.Get(i).Valid() &&
		points.Get(i).OnSurface() &&
		Dist2 (points.Get(i).P(), p0) <= xh2 &&
		invpindex.Get(i) <= 0)
	  {
		invpindex.Elem(i) =
		  locpoints.Append (points.Get(i).P());
		pindex.Append(i);
	  }
	  }
	*/

	pgeominfo.SetSize(locpoints.Size());
	for (int i = 0; i < pgeominfo.Size(); i++)
	{
	  pgeominfo[i].Init();
	}


	for (int i = 0; i < loclines.Size(); i++)
	{
	  for (int j = 0; j < 2; j++)
	  {
	  int lpi = loclines[i][j];

	  PointGeomInfo gi = lines[lindex[i]].GetGeomInfo(j + 1);
	  pgeominfo.Elem(lpi).AddPointGeomInfo(gi);

	  /*
	    if (pgeominfo.Elem(lpi).cnt == MULTIPOINTGEOMINFO_MAX)
	    break;

	    const PointGeomInfo & gi =
	    lines.Get(lindex.Get(i)).GetGeomInfo (j);

	    PointGeomInfo * pgi = pgeominfo.Elem(lpi).mgi;

	    int found = 0;
	    for (k = 0; k < pgeominfo.Elem(lpi).cnt; k++)
	    if (pgi[k].trignum == gi.trignum)
	    found = 1;

	    if (!found)
	    {
	    pgi[pgeominfo.Elem(lpi).cnt] = gi;
	    pgeominfo.Elem(lpi).cnt++;
	    }
	  */
	  }
	}

	for (int i = 0; i < locpoints.Size(); i++)
	{
	int pi = pindex[i];

	if (points[pi].mgi)
	{
	  for (int j = 1; j <= points[pi].mgi.GetNPGI(); j++)
	  {
		pgeominfo[i].AddPointGeomInfo(points[pi].mgi.GetPGI(j));
	  }
	}
	}

	if (loclines.Size() == 1)
	{
	Console.Write("loclines.Size = 1");
	Console.Write("\n");
	(*testout) << "loclines.size = 1" << "\n" << " h = " << xh << "\n" << " nearline.size = " << nearlines.Size() << "\n" << " p0 = " << p0 << "\n";
	}

	return lines[baselineindex].LineClass();
  }

  ///
  public void DeleteLine(int li)
  {
	int pi;

	nfl--;

	for (int i = 1; i <= 2; i++)
	{
	pi = lines[li].L().I(i);
	points[pi].RemoveLine();

	if (!points[pi].Valid())
	{
		delpointl.Append(pi);
		if (points[pi].mgi)
		{
		cpointsearchtree.DeleteElement(pi);
		points[pi].mgi = null;
		points[pi].mgi = null;
		}

			pointsearchtree.DeleteElement(pi);
	}
	}

	if (allflines)
	{
	allflines.Set(new INDEX_2(GetGlobalIndex(lines[li].L().I1()), GetGlobalIndex(lines[li].L().I2())), 2);
	}

	lines[li].Invalidate();
	linesearchtree.DeleteElement(li);

	dellinel.Append(li);
  }

  ///

  /*
  void AdFront2 :: GetPoints (Array<Point<3> > & apoints) const
  {
    apoints.Append (points);
    // for (int i = 0; i < points.Size(); i++)
    // apoints.Append (points[i].P());
  }
  */



  public int AddPoint(Point < 3> p, PointIndex globind, MultiPointGeomInfo mgi = null, bool pointonsurface = true)
  {
	// inserts at empty position or resizes array
	int pi;

	if (delpointl.Size() != 0)
	{
	pi = delpointl.Last();
	delpointl.DeleteLast();

	points[pi] = new FrontPoint2(p, globind, mgi, pointonsurface);
	}
	else
	{
	points.Append(new FrontPoint2(p, globind, mgi, pointonsurface));
		pi = points.Size() - 1;
	}

	if (mgi != null)
	{
	  cpointsearchtree.Insert(p, pi);
	}

	if (pointonsurface)
	{
	  pointsearchtree.Insert(p, pi);
	}

	return pi;
  }

  ///
  public int AddLine(int pi1, int pi2, PointGeomInfo gi1, PointGeomInfo gi2)
  {
	int minfn;
	int li;

	FrontPoint2 p1 = points[pi1];
	FrontPoint2 p2 = points[pi2];


	nfl++;

	p1.AddLine();
	p2.AddLine();

	minfn = netgen.GlobalMembers.min2(p1.FrontNr(), p2.FrontNr());
	p1.DecFrontNr(minfn + 1);
	p2.DecFrontNr(minfn + 1);

	if (dellinel.Size() != 0)
	{
	li = dellinel.Last();
	dellinel.DeleteLast();
	lines[li] = new FrontLine(new INDEX_2(pi1, pi2));
	}
	else
	{
	lines.Append(new FrontLine(new INDEX_2(pi1, pi2)));
		li = lines.Size() - 1;
	}


	if (!gi1.trignum || !gi2.trignum)
	{
	Console.Write("WARNING: in AdFront::AddLine, illegal geominfo");
	Console.Write("\n");
	}

	lines[li].SetGeomInfo(gi1, gi2);

	Box3d lbox = new Box3d();
	lbox.SetPoint(p1.P());
	lbox.AddPoint(p2.P());

	linesearchtree.Insert(lbox.PMin(), lbox.PMax(), li);

	if (allflines)
	{
	if (allflines.Used(new INDEX_2(GetGlobalIndex(pi1), GetGlobalIndex(pi2))))
	{
		cerr << "ERROR Adfront2::AddLine: line exists" << "\n";
		(*testout) << "ERROR Adfront2::AddLine: line exists" << "\n";
	}

	allflines.Set(new INDEX_2(GetGlobalIndex(pi1), GetGlobalIndex(pi2)), 1);
	}

	return li;
  }

  ///
  public int ExistsLine(int pi1, int pi2)
  {
	if (!allflines)
	{
	  return 0;
	}
	if (allflines.Used(new INDEX_2(pi1, pi2)))
	{
	  return allflines.Get(new INDEX_2(pi1, pi2));
	}
	else
	{
	  return 0;
	}
  }

  ///
  public void IncrementClass(int li)
  {
	lines[li].IncrementClass();
  }

  ///
  public void ResetClass(int li)
  {
	lines[li].ResetClass();
  }

  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const PointGeomInfo & GetLineGeomInfo(int li, int lend) const
  public PointGeomInfo GetLineGeomInfo(int li, int lend)
  {
		return lines[li].GetGeomInfo(lend);
  }
  ///

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: PointIndex GetGlobalIndex(int pi) const
  public PointIndex GetGlobalIndex(int pi)
  {
	return points[pi].GlobalIndex();
  }


  /// is Point p inside Surface (flat geometry only)
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool Inside(const Point<2> & p) const
  public bool Inside(Point < 2> p)
  {
	int cnt;
	Vec < 2> n;
	Vec < 3> v1;
	DenseMatrix a = new DenseMatrix(2);
	DenseMatrix ainv = new DenseMatrix(2);
	Vector b = new Vector(2);
	Vector u = new Vector(2);

	// quasi-random numbers:
	n(0) = 0.123871;
	n(1) = 0.15432;

	cnt = 0;
	for (int i = 0; i < lines.Size(); i++)
	{
	  if (lines[i].Valid())
	  {
	  const Point < 3> & p1 = points[lines[i].L().I1()].P();
	  const Point < 3> & p2 = points[lines[i].L().I2()].P();

	  v1 = p2 - p1;

	  a(0, 0) = v1(0);
	  a(1, 0) = v1(1);

	  a(0, 1) = -n(0);
	  a(1, 1) = -n(1);

	  b(0) = p(0) - p1(0);
	  b(1) = p(1) - p1(1);

	  netgen.GlobalMembers.CalcInverse(a, ainv);
	  ainv.Mult(b, u);

	  if (u(0) >= 0 && u(0) <= 1 && u(1) > 0)
	  {
		cnt++;
	  }
	  }
	}

	return ((cnt % 2) != 0);
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool SameSide(const Point<2> & lp1, const Point<2> & lp2, const Array<int> * testfaces = null) const
  public bool SameSide(Point < 2> lp1, Point < 2> lp2, Array<int> testfaces = null)
  {
	int cnt = 0;

	if (testfaces != null)
	{
		for (int ii = 0; ii < testfaces.Size(); ii++)
		{
		  if (lines[testfaces[ii]].Valid())
		  {
			  int i = testfaces[ii];
			  const Point < 3> & p13d = points[lines[i].L().I1()].P();
			  const Point < 3> & p23d = points[lines[i].L().I2()].P();

			  Point < 2> p1(p13d(0), p13d(1));
			  Point < 2> p2(p23d(0), p23d(1));

			  // p1 + alpha v = lp1 + beta vl
			  Vec < 2> v = p2 - p1;
			  Vec < 2> vl = lp2 - lp1;
			  Mat < 2,2> mat, inv;
			  Vec < 2> rhs, sol;
			  mat(0,0) = v(0);
			  mat(1,0) = v(1);
			  mat(0,1) = -vl(0);
			  mat(1,1) = -vl(1);
			  rhs = lp1 - p1;

			  if (Det(mat) == 0)
			  {
				  continue;
			  }
			  netgen.GlobalMembers.CalcInverse(mat, inv);
			  sol = inv * rhs;
			  if ((sol(0) >= 0) && (sol(0) <= 1) && (sol(1) >= 0) && (sol(1) <= 1))
			  {
					cnt++;
			  }
		  }
		}

	}
	else
	{
		for (int i = 0; i < lines.Size(); i++)
		{
		  if (lines[i].Valid())
		  {
			  const Point < 3> & p13d = points[lines[i].L().I1()].P();
			  const Point < 3> & p23d = points[lines[i].L().I2()].P();

			  Point < 2> p1(p13d(0), p13d(1));
			  Point < 2> p2(p23d(0), p23d(1));

			  // p1 + alpha v = lp1 + beta vl
			  Vec < 2> v = p2 - p1;
			  Vec < 2> vl = lp2 - lp1;
			  Mat < 2,2> mat, inv;
			  Vec < 2> rhs, sol;
			  mat(0,0) = v(0);
			  mat(1,0) = v(1);
			  mat(0,1) = -vl(0);
			  mat(1,1) = -vl(1);
			  rhs = lp1 - p1;

			  if (Det(mat) == 0)
			  {
				  continue;
			  }
			  netgen.GlobalMembers.CalcInverse(mat, inv);
			  sol = inv * rhs;
			  if ((sol(0) >= 0) && (sol(0) <= 1) && (sol(1) >= 0) && (sol(1) <= 1))
			  {
					cnt++;
			  }
		  }
		}
	}
	return ((cnt % 2) == 0);
  }

  /*
  {
    return Inside (lp1) == Inside (lp2);
  }
  */

  ///
  public void SetStartFront()
  {
	for (int i = lines.Begin(); i < lines.End(); i++)
	{
	  if (lines[i].Valid())
	  {
	for (int j = 1; j <= 2; j++)
	{
	  points[lines[i].L().I(j)].DecFrontNr(0);
	}
	  }
	}
  }

  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void PrintOpenSegments(ostream & ost) const
  public void PrintOpenSegments(ostream ost)
  {
	if (nfl > 0)
	{
	ost << nfl << " open front segments left:" << "\n";
	for (int i = lines.Begin(); i < lines.End(); i++)
	{
	  if (lines[i].Valid())
	  {
		ost << i << ": " << GetGlobalIndex(lines[i].L().I1()) << "-" << GetGlobalIndex(lines[i].L().I2()) << "\n";
	  }
	}
	}
  }
}









