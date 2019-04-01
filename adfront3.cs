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
/* File:   adfront3.hh                                                    */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Okt. 95                                                    */
/**************************************************************************/

/*
  Advancing front class for volume meshing
*/



/// Point in advancing front
public class FrontPoint3
{
  /// coordinates
  private Point < 3> p;
  /// global node index
  private PointIndex globalindex = new PointIndex();
  /// number of faces connected to point 
  private int nfacetopoint;
  /// distance to original boundary
  private int frontnr;
  /// 
  private int cluster;
  ///
  public FrontPoint3()
  {
	globalindex.Invalidate(); //  = -1;
	nfacetopoint = 0;
	frontnr = 1000;
	cluster = 0;
  }

  ///
  public FrontPoint3(Point < 3> ap, PointIndex agi)
  {
	p = ap;
	globalindex = agi;
	nfacetopoint = 0;
	frontnr = 1000;
	cluster = 0;
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
//ORIGINAL LINE: PointIndex GlobalIndex() const
  public PointIndex GlobalIndex()
  {
	  return globalindex;
  }

  ///
  public void AddFace()
  {
	  nfacetopoint++;
  }

  /// if last face is removed, then point is invalidated
  public void RemoveFace()
  {
	nfacetopoint--;
	if (nfacetopoint == 0)
	{
		nfacetopoint = -1;
	}
  }

  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool Valid() const
  public bool Valid()
  {
	  return nfacetopoint >= 0;
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

  ///
//C++ TO C# CONVERTER TODO TASK: C# has no concept of a 'friend' class:
//  friend class AdFront3;
}



public class MiniElement2d
{
  protected int np;
  protected PointIndex[] pnum = Arrays.InitializeWithDefaultInstances<PointIndex>(4); // can be global or local nums
  protected bool deleted;
  public MiniElement2d()
  {
	  np = 3;
	  deleted = 0;
  }
  public MiniElement2d(int anp)
  {
	  np = anp;
	  deleted = 0;
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetNP() const
  public int GetNP()
  {
	  return np;
  }
  public PointIndex this [int i]
  {
	  get
	  {
		  return pnum[i];
	  }
	  set
	  {
		  pnum[i] = value;
	  }
  }
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const PointIndex operator [] (int i) const
  public PointIndex this [int i]
  {
	  get
	  {
		  return pnum[i];
	  }
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const PointIndex PNum(int i) const
  public PointIndex PNum(int i)
  {
	  return pnum[i - 1];
  }
  public PointIndex PNum(int i)
  {
	  return pnum[i - 1];
  }
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const PointIndex PNumMod(int i) const
  public PointIndex PNumMod(int i)
  {
	  return pnum[(i - 1) % np];
  }
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: auto PNums() const
//C++ TO C# CONVERTER TODO TASK: The return type of the following function could not be determined:
  public auto PNums()
  {
	  return new FlatArray<const PointIndex> (np, pnum[0]);
  }
  public void Delete()
  {
	  deleted = true;
	  foreach (PointIndex p in pnum)
	  {
		  p.Invalidate();
	  }
  }
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool IsDeleted() const
  public bool IsDeleted()
  {
	  return deleted;
  }
}




/// Face in advancing front
public class FrontFace
{
  ///
  private MiniElement2d f = new MiniElement2d();
  ///
  private int qualclass;
  ///
  private char oldfront;
  ///
  private int hashvalue;
  ///
  private int cluster;

  ///

  /* ********************** FrontFace ********************** */

  public FrontFace()
  {
	qualclass = 1;
	oldfront = 0;
	hashvalue = 0;
	cluster = 0;
  }

  ///
  public FrontFace(MiniElement2d af)
  {
	f = af;
	oldfront = 0;
	qualclass = 1;
	hashvalue = 0;
  }

  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const MiniElement2d & Face() const
  public MiniElement2d Face()
  {
	  return f;
  }

  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int QualClass() const
  public int QualClass()
  {
	  return qualclass;
  }

  ///
  public void IncrementQualClass()
  {
	  qualclass++;
  }

  ///
  public void ResetQualClass()
  {
	if (qualclass > 1)
	{
	qualclass = 1;
	oldfront = 0;
	}
  }

  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool Valid() const
  public bool Valid()
  {
	  return !f.IsDeleted();
  }

  ///
  public void Invalidate()
  {
	f.Delete();
	oldfront = 0;
	qualclass = 1000;
  }

  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int HashValue() const
  public int HashValue()
  {
	  return hashvalue;
  }

  ///
  public void SetHashValue(int hv)
  {
	  hashvalue = hv;
  }

  ///
//C++ TO C# CONVERTER TODO TASK: C# has no concept of a 'friend' class:
//  friend class AdFront3;

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int Cluster() const
  public int Cluster()
  {
	  return cluster;
  }
}




/// Advancing front, 3D.
public class AdFront3
{
  ///
  private Array<FrontPoint3, PointIndex.BASE, PointIndex> points = new Array<FrontPoint3, PointIndex.BASE, PointIndex>();
  ///
  private Array<FrontFace> faces = new Array<FrontFace>();
  ///
  private Array<PointIndex> delpointl = new Array<PointIndex>();

  /// which points are connected to pi ?
  private TABLE<int, PointIndex.BASE> connectedpairs;

  /// number of total front faces;
  private int nff;
  /// number of quads in front
  private int nff4;

  ///
  private double vol;

  ///
  private GeomSearch3d hashtable = new GeomSearch3d();

  /// 
  private int hashon;

  ///
  private int hashcreated;

  /// counter for rebuilding internal tables
  private int rebuildcounter;
  /// last base element
  private int lasti;
  /// minimal selection-value of baseelements
  private int minval;
  private Array<PointIndex, PointIndex.BASE, PointIndex> invpindex = new Array<PointIndex, PointIndex.BASE, PointIndex>();
  private Array<char> pingroup = new Array<char>();

  ///
  private BoxTree < 3> * facetree;

  ///

  /* ********************** AddFront ********************** */


  public AdFront3()
  {
	nff = 0;
	nff4 = 0;
	vol = 0;

	hashon = 1;
	hashcreated = 0;
	if (hashon)
	{
	  hashtable.Init(points, faces);
	}

	facetree = null;
	connectedpairs = null;

	rebuildcounter = -1;
	lasti = 0;
	minval = -1;
  }

  ///
  public void Dispose()
  {
	if (facetree != null)
	{
		facetree.Dispose();
	}
	if (connectedpairs != null)
	{
		connectedpairs.Dispose();
	}
  }

  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetPoints(Array<Point<3>> & apoints) const
  public void GetPoints(Array<Point < 3>> apoints)
  {
	for (PointIndex pi = points.Begin(); pi < points.End(); pi++)
	{

	  apoints.Append(points[pi].P());
	}
  }

  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetNP() const
  public int GetNP()
  {
	  return points.Size();
  }
  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const Point<3> & GetPoint(PointIndex pi) const
  public Point < 3> GetPoint(PointIndex pi)
  {
	  return points[pi].P();
  }
  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetNF() const
  public int GetNF()
  {
	  return nff;
  }
  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const MiniElement2d & GetFace(int i) const
  public MiniElement2d GetFace(int i)
  {
	  return faces.Get(i).Face();
  }
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const auto & Faces() const
//C++ TO C# CONVERTER TODO TASK: The return type of the following function could not be determined:
  public auto Faces()
  {
	  return faces;
  }
  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void Print() const;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  void Print();
  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool Empty() const
  public bool Empty()
  {
	  return nff == 0;
  }
  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool Empty(int elnp) const
  public bool Empty(int elnp)
  {
	if (elnp == 4)
	{
	  return (nff4 == 0);
	}
	return (nff - nff4 == 0);
  }
  ///
  public int SelectBaseElement()
  {
	int i;
	int hi;
	int fstind;

	/*
	static int minval = -1;
	static int lasti = 0;
	static int counter = 0;
	*/
	if (rebuildcounter <= 0)
	{
		RebuildInternalTables();
		rebuildcounter = nff / 10 + 1;

		lasti = 0;
	}
	rebuildcounter--;

	/*
	if (faces.Size() > 2 * nff)
	  {
	    // compress facelist
  
	    RebuildInternalTables ();
	    lasti = 0;
	  }
	  */

	fstind = 0;

	for (i = lasti + 1; i <= faces.Size() && !fstind; i++)
	{
	  if (faces.Elem(i).Valid())
	  {
	  hi = faces.Get(i).QualClass() + points[faces.Get(i).Face().PNum(1)].FrontNr() + points[faces.Get(i).Face().PNum(2)].FrontNr() + points[faces.Get(i).Face().PNum(3)].FrontNr();

	  if (hi <= minval)
	  {
		  minval = hi;
		  fstind = i;
		  lasti = fstind;
	  }
	  }
	}

	if (fstind == 0)
	{
		minval = INT_MAX;
		for (i = 1; i <= faces.Size(); i++)
		{
	  if (faces.Elem(i).Valid())
	  {
		  hi = faces.Get(i).QualClass() + points[faces.Get(i).Face().PNum(1)].FrontNr() + points[faces.Get(i).Face().PNum(2)].FrontNr() + points[faces.Get(i).Face().PNum(3)].FrontNr();

		  if (hi <= minval)
		  {
		  minval = hi;
		  fstind = i;
		  lasti = 0;
		  }
	  }
		}
	}


	return fstind;
  }

  ///
  public void CreateTrees()
  {
	int i;
	int j;
	PointIndex pi = new PointIndex();
	Point3d pmin = new Point3d();
	Point3d pmax = new Point3d();

	for (pi = PointIndex.BASE; pi < GetNP() + PointIndex.BASE; pi++)
	{
		const Point < 3> & p = GetPoint(pi);
		if (pi == PointIndex.BASE)
		{
		pmin = p;
		pmax = p;
		}
		else
		{
		pmin.SetToMin(p);
		pmax.SetToMax(p);
		}
	}

//C++ TO C# CONVERTER TODO TASK: The following line was determined to be a copy assignment (rather than a reference assignment) - this should be verified and a 'CopyFrom' method should be created:
//ORIGINAL LINE: pmax = pmax + 0.5 * (pmax - pmin);
	pmax.CopyFrom(pmax + 0.5 * (pmax - pmin));
//C++ TO C# CONVERTER TODO TASK: The following line was determined to be a copy assignment (rather than a reference assignment) - this should be verified and a 'CopyFrom' method should be created:
//ORIGINAL LINE: pmin = pmin + 0.5 * (pmin - pmax);
	pmin.CopyFrom(pmin + 0.5 * (pmin - pmax));

	facetree = null;
	facetree = new BoxTree < 3> (pmin, pmax);

	for (i = 1; i <= GetNF(); i++)
	{
		MiniElement2d el = GetFace(i);
		pmin = GetPoint(el[0]);
//C++ TO C# CONVERTER TODO TASK: The following line was determined to be a copy assignment (rather than a reference assignment) - this should be verified and a 'CopyFrom' method should be created:
//ORIGINAL LINE: pmax = pmin;
		pmax.CopyFrom(pmin);
		for (j = 1; j < 3; j++)
		{
		const Point < 3> & p = GetPoint(el[j]);
		pmin.SetToMin(p);
		pmax.SetToMax(p);
		}
//C++ TO C# CONVERTER TODO TASK: The following line was determined to be a copy assignment (rather than a reference assignment) - this should be verified and a 'CopyFrom' method should be created:
//ORIGINAL LINE: pmax = pmax + 0.01 * (pmax - pmin);
		pmax.CopyFrom(pmax + 0.01 * (pmax - pmin));
//C++ TO C# CONVERTER TODO TASK: The following line was determined to be a copy assignment (rather than a reference assignment) - this should be verified and a 'CopyFrom' method should be created:
//ORIGINAL LINE: pmin = pmin + 0.01 * (pmin - pmax);
		pmin.CopyFrom(pmin + 0.01 * (pmin - pmax));
		//      (*testout) << "insert " << i << ": " << pmin << " - " << pmax << "\n";
		facetree.Insert(pmin, pmax, i);
	}
  }

  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetIntersectingFaces(const Point<3> & pmin, const Point<3> & pmax, Array<int> & ifaces) const
  public void GetIntersectingFaces(Point < 3> pmin, Point < 3> pmax, Array<int> ifaces)
  {
	facetree.GetIntersecting(pmin, pmax, ifaces);
  }

  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetFaceBoundingBox(int i, Box3d & box) const
  public void GetFaceBoundingBox(int i, Box3d box)
  {
	FrontFace face = faces.Get(i);
	box.SetPoint(points[face.f[0]].p);
	box.AddPoint(points[face.f[1]].p);
	box.AddPoint(points[face.f[2]].p);
  }

  ///
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  int GetLocals(int baseelement, Array<Point3d, PointIndex::BASE> locpoints, Array<MiniElement2d> locfaces, Array<PointIndex, PointIndex::BASE> pindex, Array<INDEX> findex, INDEX_2_HASHTABLE<int> connectedpairs, float xh, float relh, INDEX facesplit);

  ///

  // returns all points connected with fi
  public void GetGroup(int fi, Array<MeshPoint, PointIndex.BASE> grouppoints, Array<MiniElement2d> groupelements, Array<PointIndex, PointIndex.BASE> pindex, Array<int> findex)
  {
	// static Array<char> pingroup;
	int changed;

	pingroup.SetSize(points.Size());

	pingroup = 0;
	for (int j = 1; j <= 3; j++)
	{
	  pingroup.Elem(faces.Get(fi).Face().PNum(j)) = 1;
	}

	do
	{
		changed = 0;

		/*
		for (i = 1; i <= faces.Size(); i++)
	  if (faces.Get(i).Valid())
		{
		  const MiniElement2d & face = faces.Get(i).Face();
  
		  int fused = 0;
		  for (j = 1; j <= 3; j++)
			if (pingroup.Elem(face.PNum(j)))
		  fused++;
		
		  if (fused >= 2)
			for (j = 1; j <= 3; j++)
		  if (!pingroup.Elem(face.PNum(j)))
			{
			  pingroup.Elem(face.PNum(j)) = 1;
			  changed = 1;
			}
		}
		*/
		foreach (var f in faces)
		{
	  if (f.Valid())
	  {
		  MiniElement2d face = f.Face();

		  int fused = 0;
		  for (int j = 1; j <= 3; j++)
		  {
			if (pingroup.Elem(face.PNum(j)))
			{
		  fused++;
			}
		  }

		  if (fused >= 2)
		  {
			for (int j = 1; j <= 3; j++)
			{
		  if (!pingroup.Elem(face.PNum(j)))
		  {
			  pingroup.Elem(face.PNum(j)) = 1;
			  changed = 1;
		  }
			}
		  }
	  }
		}

	} while (changed != 0);

	invpindex.SetSize(points.Size());

	// for (PointIndex pi = points.Begin(); pi < points.End(); pi++)
	foreach (PointIndex pi in points.Range())
	{
	  if (points[pi].Valid())
	  {
	  grouppoints.Append(points[pi].P());
	  pindex.Append(pi);
	  invpindex[pi] = pindex.Size();
	  }
	}

	for (int i = 1; i <= faces.Size(); i++)
	{
	  if (faces.Get(i).Valid())
	  {
	  int fused = 0;
	  for (int j = 1; j <= 3; j++)
	  {
		if (pingroup.Get(faces.Get(i).Face().PNum(j)))
		{
		  fused++;
		}
	  }

	  if (fused >= 2)
	  {
		  groupelements.Append(faces.Get(i).Face());
		  findex.Append(i);
	  }
	  }
	}

	/*
	for (int i = 1; i <= groupelements.Size(); i++)
	  for (int j = 1; j <= 3; j++)
	    {
	  groupelements.Elem(i).PNum(j) =
		invpindex.Get(groupelements.Elem(i).PNum(j));
	    }
	*/
	foreach (var e in groupelements)
	{
	  for (int j = 1; j <= 3; j++)
	  {
		e.PNum(j) = invpindex.Get(e.PNum(j));
	  }
	}
  }

  ///
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  void DeleteFace(INDEX fi);
  ///
  public PointIndex AddPoint(Point < 3> p, PointIndex globind)
  {
	if (delpointl.Size())
	{
		PointIndex pi = delpointl.Last();
		delpointl.DeleteLast();

		points[pi] = new FrontPoint3(p, globind);
		return new PointIndex(pi);
	}
	else
	{
		points.Append(new FrontPoint3(p, globind));
		return --points.End();
		// return points.Size()-1+PointIndex::BASE;
	}
  }

  ///
  public int AddFace(MiniElement2d aface)
  {
	int i;
	int minfn;

	nff++;

	for (i = 0; i < aface.GetNP(); i++)
	{
	  points[aface[i]].AddFace();
	}

	Point3d p1 = points[aface[0]].P();
	Point3d p2 = points[aface[1]].P();
	Point3d p3 = points[aface[2]].P();

	vol += 1.0 / 6.0 * (p1.X() + p2.X() + p3.X()) * ((p2.Y() - p1.Y()) * (p3.Z() - p1.Z()) - (p2.Z() - p1.Z()) * (p3.Y() - p1.Y()));

	if (aface.GetNP() == 4)
	{
		nff4++;
		Point3d p4 = points[aface[3]].P();
		vol += 1.0 / 6.0 * (p1.X() + p3.X() + p4.X()) * ((p3.Y() - p1.Y()) * (p4.Z() - p1.Z()) - (p3.Z() - p1.Z()) * (p4.Y() - p1.Y()));
	}


	minfn = 1000;
	for (i = 0; i < aface.GetNP(); i++)
	{
		int fpn = points[aface[i]].FrontNr();
		if (i == 0 || fpn < minfn)
		{
	  minfn = fpn;
		}
	}


	int cluster = 0;
	for (i = 1; i <= aface.GetNP(); i++)
	{
		if (points[aface.PNum(i)].cluster)
		{
	  cluster = points[aface.PNum(i)].cluster;
		}
	}
	for (i = 1; i <= aface.GetNP(); i++)
	{
	  points[aface.PNum(i)].cluster = cluster;
	}


	for (i = 1; i <= aface.GetNP(); i++)
	{
	  points[aface.PNum(i)].DecFrontNr(minfn + 1);
	}

	faces.Append(new FrontFace(aface));
	int nfn = faces.Size();
	faces.Elem(nfn).cluster = cluster;

	if (hashon && hashcreated)
	{
	  hashtable.AddElem(aface, nfn);
	}

	return nfn;
  }

  ///
  public int AddConnectedPair(INDEX_2 apair)
  {
	if (!connectedpairs)
	{
	  connectedpairs = new TABLE<int, PointIndex.BASE> (GetNP());
	}

	connectedpairs.Add(apair.I1(), apair.I2());
	connectedpairs.Add(apair.I2(), apair.I1());

	return 0;
  }

  ///
  public void IncrementClass(INDEX fi)
  {
	  faces.Elem(fi).IncrementQualClass();
  }

  ///
  public void ResetClass(INDEX fi)
  {
	  faces.Elem(fi).ResetQualClass();
  }

  ///
  public void SetStartFront(int baseelnp = 0)
  {
	for (int i = 1; i <= faces.Size(); i++)
	{
	  if (faces.Get(i).Valid())
	  {
	  MiniElement2d face = faces.Get(i).Face();
	  for (int j = 1; j <= 3; j++)
	  {
		points[face.PNum(j)].DecFrontNr(0);
	  }
	  }
	}

	/*
	if (baseelnp)
	  {
	    for (i = 1; i <= faces.Size(); i++)
	  if (faces.Get(i).Valid() && faces.Get(i).Face().GetNP() != baseelnp)
		faces.Elem(i).qualclass = 1000;
	  }
	  */
  }

  /// is Point p inside Surface ?
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool Inside(const Point<3> & p) const
  public bool Inside(Point < 3> p)
  {
	int cnt;
	Vec3d n = new Vec3d();
	Vec3d v1 = new Vec3d();
	Vec3d v2 = new Vec3d();
	DenseMatrix a = new DenseMatrix(3);
	DenseMatrix ainv = new DenseMatrix(3);
	Vector b = new Vector(3);
	Vector u = new Vector(3);

	// random numbers:
	n.X() = 0.123871;
	n.Y() = 0.15432;
	n.Z() = -0.43989;

	cnt = 0;
	for (int i = 1; i <= faces.Size(); i++)
	{
	  if (faces.Get(i).Valid())
	  {
	  const Point < 3> & p1 = points[faces.Get(i).Face().PNum(1)].P();
	  const Point < 3> & p2 = points[faces.Get(i).Face().PNum(2)].P();
	  const Point < 3> & p3 = points[faces.Get(i).Face().PNum(3)].P();

	  v1 = p2 - p1;
	  v2 = p3 - p1;

	  a(0, 0) = v1.X();
	  a(1, 0) = v1.Y();
	  a(2, 0) = v1.Z();
	  a(0, 1) = v2.X();
	  a(1, 1) = v2.Y();
	  a(2, 1) = v2.Z();
	  a(0, 2) = -n.X();
	  a(1, 2) = -n.Y();
	  a(2, 2) = -n.Z();

	  b(0) = p(0) - p1(0);
	  b(1) = p(1) - p1(1);
	  b(2) = p(2) - p1(2);

	  netgen.GlobalMembers.CalcInverse(a, ainv);
	  ainv.Mult(b, u);

	  if (u(0) >= 0 && u(1) >= 0 && u(0) + u(1) <= 1 && u(2) > 0)
	  {
		  cnt++;
	  }
	  }
	}

	return ((cnt % 2) != 0);
  }

  /// both points on same side ?
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int SameSide(const Point<3> & lp1, const Point<3> & lp2, const Array<int> * testfaces = null) const
  public int SameSide(Point < 3> lp1, Point < 3> lp2, Array<int> testfaces = null)
  {
	const Point < 3> *line[2];
	line[0] = lp1;
	line[1] = lp2;


	Point3d pmin = new Point3d(lp1);
	Point3d pmax = new Point3d(lp1);
	pmin.SetToMin(lp2);
	pmax.SetToMax(lp2);

	ArrayMem<int, 100> aprif = new ArrayMem<int, 100>();
	aprif.SetSize(0);

	if (testfaces == null)
	{
	  facetree.GetIntersecting(pmin, pmax, aprif);
	}
	else
	{
	  for (int i = 1; i <= testfaces.Size(); i++)
	  {
		aprif.Append(testfaces.Get(i));
	  }
	}

	int cnt = 0;
	for (int ii = 1; ii <= aprif.Size(); ii++)
	{
		int i = aprif.Get(ii);

		if (faces.Get(i).Valid())
		{
		const Point < 3> *tri[3];
		tri[0] = &points[faces.Get(i).Face().PNum(1)].P();
		tri[1] = &points[faces.Get(i).Face().PNum(2)].P();
		tri[2] = &points[faces.Get(i).Face().PNum(3)].P();

		if (IntersectTriangleLine(tri[0], line[0]))
		{
		  cnt++;
		}
		}
	}

	return ((cnt + 1) % 2);
  }


  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: PointIndex GetGlobalIndex(PointIndex pi) const
  public PointIndex GetGlobalIndex(PointIndex pi)
  {
	  return points[pi].GlobalIndex();
  }
  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: double Volume() const
  public double Volume()
  {
	  return vol;
  }


//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private int RebuildInternalTables_timer_a = NgProfiler.CreateTimer("Adfront3::RebuildInternal A");
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private int RebuildInternalTables_timer_b = NgProfiler.CreateTimer("Adfront3::RebuildInternal B");
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private int RebuildInternalTables_timer_c = NgProfiler.CreateTimer("Adfront3::RebuildInternal C");
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private int RebuildInternalTables_timer_d = NgProfiler.CreateTimer("Adfront3::RebuildInternal D");

  private void RebuildInternalTables()
  {
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static int timer_a = NgProfiler::CreateTimer("Adfront3::RebuildInternal A");
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static int timer_b = NgProfiler::CreateTimer("Adfront3::RebuildInternal B");
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static int timer_c = NgProfiler::CreateTimer("Adfront3::RebuildInternal C");
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static int timer_d = NgProfiler::CreateTimer("Adfront3::RebuildInternal D");


	NgProfiler.StartTimer(RebuildInternalTables_timer_a);
	int hi = 0;
	for (int i = 1; i <= faces.Size(); i++)
	{
	  if (faces.Get(i).Valid())
	  {
	  hi++;
	  if (hi < i)
	  {
		faces.Elem(hi) = faces.Get(i);
	  }
	  }
	}

	faces.SetSize(nff);

	int np = points.Size();

	for (PointIndex pi = points.Begin(); pi < points.End(); pi++)
	{
	  points[pi].cluster = pi;
	}

	NgProfiler.StopTimer(RebuildInternalTables_timer_a);
	NgProfiler.StartTimer(RebuildInternalTables_timer_b);

	int change;
	do
	{
		change = 0;
		for (int i = 1; i <= faces.Size(); i++)
		{
		MiniElement2d el = faces.Get(i).Face();

		int mini = points[el.PNum(1)].cluster;
		int maxi = mini;

		for (int j = 2; j <= 3; j++)
		{
			int ci = points[el.PNum(j)].cluster;
			if (ci < mini)
			{
				mini = ci;
			}
			if (ci > maxi)
			{
				maxi = ci;
			}
		}

		if (mini < maxi)
		{
			change = 1;
			for (int j = 1; j <= 3; j++)
			{
		  points[el.PNum(j)].cluster = mini;
			}
		}
		}
	} while (change != 0);


	NgProfiler.StopTimer(RebuildInternalTables_timer_b);
	NgProfiler.StartTimer(RebuildInternalTables_timer_c);




	BitArrayChar<PointIndex.BASE> usecl = new BitArrayChar<PointIndex.BASE>(np);
	usecl.Clear();
	for (int i = 1; i <= faces.Size(); i++)
	{
		usecl.Set(points[faces.Get(i).Face().PNum(1)].cluster);
		faces.Elem(i).cluster = points[faces.Get(i).Face().PNum(1)].cluster;
	}
	int cntcl = 0;
	for (int i = PointIndex.BASE; i < np + PointIndex.BASE; i++)
	{
	  if (usecl.Test(i))
	  {
		cntcl++;
	  }
	}

	Array<double, PointIndex.BASE> clvol = new Array<double, PointIndex.BASE>(np);
	clvol = 0.0;

	for (int i = 1; i <= faces.Size(); i++)
	{
		MiniElement2d face = faces.Get(i).Face();

		Point3d p1 = points[face.PNum(1)].P();
		Point3d p2 = points[face.PNum(2)].P();
		Point3d p3 = points[face.PNum(3)].P();

		double vi = 1.0 / 6.0 * (p1.X() + p2.X() + p3.X()) * ((p2.Y() - p1.Y()) * (p3.Z() - p1.Z()) - (p2.Z() - p1.Z()) * (p3.Y() - p1.Y()));

		if (face.GetNP() == 4)
		{
		Point3d p4 = points[face.PNum(4)].P();
		vi += 1.0 / 6.0 * (p1.X() + p3.X() + p4.X()) * ((p3.Y() - p1.Y()) * (p4.Z() - p1.Z()) - (p3.Z() - p1.Z()) * (p4.Y() - p1.Y()));
		}

		clvol[faces.Get(i).cluster] += vi;
	}

	NgProfiler.StopTimer(RebuildInternalTables_timer_c);
	NgProfiler.StartTimer(RebuildInternalTables_timer_d);



	int negvol = 0;
	for (int i = PointIndex.BASE; i < clvol.Size() + PointIndex.BASE; i++)
	{
		if (clvol[i] < 0)
		{
	  negvol = 1;
		}
	}

	if (negvol != 0)
	{
		for (int i = 1; i <= faces.Size(); i++)
		{
	  faces.Elem(i).cluster = 1;
		}
		for (PointIndex pi = points.Begin(); pi < points.End(); pi++)
		{
	  points[pi].cluster = 1;
		}
	}

	if (hashon)
	{
	  hashtable.Create();
	}

	NgProfiler.StopTimer(RebuildInternalTables_timer_d);
  }
}







/* ********************** FrontPoint ********************** */

