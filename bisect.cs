using System;

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

public class BisectionOptions
{
  public readonly string outfilename;
  public readonly string mlfilename;
  public readonly string refinementfilename;
  public readonly string femcode;
  public int maxlevel;
  public int usemarkedelements;
  public bool refine_hp;
  public bool refine_p;
  public TaskManager task_manager = &DummyTaskManager;
  public Tracer tracer = &DummyTracer;
  public BisectionOptions()
  {
	outfilename = null;
	mlfilename = null;
	refinementfilename = null;
	femcode = null;
	maxlevel = 50;
	usemarkedelements = 0;
	refine_hp = 0;
	refine_p = 0;
  }
}

public class ZRefinementOptions
{
  public int minref;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  DLL_HEADER ZRefinementOptions();
}





//C++ TO C# CONVERTER WARNING: The original type declaration contained unconverted modifiers:
//ORIGINAL LINE: class DLL_HEADER Refinement
public class Refinement
{
  private MeshOptimize2d optimizer2d;

  public Refinement()
  {
	optimizer2d = null;
  }

  public virtual void Dispose()
  {
	;
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void Refine(Mesh & mesh) const;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  void Refine(Mesh mesh);
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  void Refine(Mesh mesh);
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void Bisect(Mesh & mesh, class BisectionOptions & opt, Array<double> * quality_loss = null) const;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  void Bisect(Mesh mesh, BisectionOptions opt, Array<double> quality_loss = null);

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void MakeSecondOrder(Mesh & mesh) const;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  void MakeSecondOrder(Mesh mesh);
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  void MakeSecondOrder(Mesh mesh);

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual void PointBetween(const Point<3> & p1, const Point<3> & p2, double secpoint, int surfi, const PointGeomInfo & gi1, const PointGeomInfo & gi2, Point<3> & newp, PointGeomInfo & newgi) const
  public virtual void PointBetween(Point < 3> p1, Point < 3> p2, double secpoint, int surfi, PointGeomInfo gi1, PointGeomInfo gi2, ref Point < 3> newp, PointGeomInfo newgi)
  {
	newp = p1 + secpoint * (p2 - p1);
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual void PointBetween(const Point<3> & p1, const Point<3> & p2, double secpoint, int surfi1, int surfi2, const EdgePointGeomInfo & ap1, const EdgePointGeomInfo & ap2, Point<3> & newp, EdgePointGeomInfo & newgi) const
  public virtual void PointBetween(Point < 3> p1, Point < 3> p2, double secpoint, int surfi1, int surfi2, EdgePointGeomInfo ap1, EdgePointGeomInfo ap2, ref Point < 3> newp, EdgePointGeomInfo newgi)
  {
	//cout << "base class edge point between" << endl;
	newp = p1 + secpoint * (p2 - p1);
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual Vec<3> GetTangent(const Point<3> & p, int surfi1, int surfi2, const EdgePointGeomInfo & ap1) const
  public virtual Vec < 3> GetTangent(Point < 3> p, int surfi1, int surfi2, EdgePointGeomInfo ap1)
  {
	cerr << "Refinement::GetTangent not overloaded" << "\n";
	return Vec < 3> (0,0,0);
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual Vec<3> GetNormal(const Point<3> & p, int surfi1, const PointGeomInfo & gi) const
  public virtual Vec < 3> GetNormal(Point < 3> p, int surfi1, PointGeomInfo gi)
  {
	cerr << "Refinement::GetNormal not overloaded" << "\n";
	return Vec < 3> (0,0,0);
  }


//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual void ProjectToSurface(Point<3> & p, int surfi) const
  public virtual void ProjectToSurface(Point < 3> p, int surfi)
  {
	if (printmessage_importance > 0)
	{
	  cerr << "Refinement :: ProjectToSurface    ERROR: no geometry set" << "\n";
	}
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual void ProjectToSurface(Point<3> & p, int surfi, const PointGeomInfo &) const
  public virtual void ProjectToSurface(Point < 3> p, int surfi, PointGeomInfo UnnamedParameter)
  {
	ProjectToSurface(p, surfi);
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual void ProjectToEdge(Point<3> & p, int surfi1, int surfi2, const EdgePointGeomInfo & egi) const
  public virtual void ProjectToEdge(Point < 3> p, int surfi1, int surfi2, EdgePointGeomInfo egi)
  {
	cerr << "Refinement::ProjectToEdge not overloaded" << "\n";
  }


//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  void ValidateSecondOrder(Mesh mesh);
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  void ValidateRefinedMesh(Mesh mesh, Array<INDEX_2> parents);

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: MeshOptimize2d * Get2dOptimizer() const
  public MeshOptimize2d Get2dOptimizer()
  {
	return optimizer2d;
  }
  public void Set2dOptimizer(MeshOptimize2d opti)
  {
	optimizer2d = opti;
  }


//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual void LocalizeEdgePoints(Mesh &) const
  public virtual void LocalizeEdgePoints(Mesh UnnamedParameter)
  {
	  ;
  }
}





namespace netgen
{
//C++ TO C# CONVERTER NOTE: C# has no need of forward class declarations:
//  class MarkedTet;
//C++ TO C# CONVERTER NOTE: C# has no need of forward class declarations:
//  class MarkedPrism;
//C++ TO C# CONVERTER NOTE: C# has no need of forward class declarations:
//  class MarkedIdentification;
//C++ TO C# CONVERTER NOTE: C# has no need of forward class declarations:
//  class MarkedTri;
//C++ TO C# CONVERTER NOTE: C# has no need of forward class declarations:
//  class MarkedQuad;




  public class MarkedTet
  {
	/// pnums of tet
	public PointIndex[] pnums = Arrays.InitializeWithDefaultInstances<PointIndex>(4);
	/// material number
	public int matindex;
	/// element marked for refinement
	/// marked = 1: marked by element marker, marked = 2 due to closure
//C++ TO C# CONVERTER TODO TASK: C# does not allow bit fields:
	public uint marked:2;
	/// flag of Arnold-Mukherjee algorithm
//C++ TO C# CONVERTER TODO TASK: C# does not allow bit fields:
	public uint flagged:1;
	/// tetedge (local coordinates 0..3)
//C++ TO C# CONVERTER TODO TASK: C# does not allow bit fields:
	public uint tetedge1:3;
//C++ TO C# CONVERTER TODO TASK: C# does not allow bit fields:
	public uint tetedge2:3;
	// marked edge of faces
	// face_j : face without node j,
	// mark_k : edge without node k

	public string faceedges = new string(new char[4]);
	// unsigned char faceedges[4];
	public bool incorder;
//C++ TO C# CONVERTER TODO TASK: C# does not allow bit fields:
	public uint order:6;

//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//	MarkedTet() = default;
	/*
	{ 
	  for (int i = 0; i < 4; i++) { faceedges[i] = 127; }
	}
	*/
//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//	MarkedTet(const MarkedTet&) = default;
//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//	MarkedTet(MarkedTet &&) = default;
//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//	MarkedTet & operator = (const MarkedTet&) = default;
//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//	MarkedTet & operator = (MarkedTet&&) = default;

  }

  public class MarkedPrism
  {
	/// 6 point numbers
	public PointIndex[] pnums = Arrays.InitializeWithDefaultInstances<PointIndex>(6);
	/// material number
	public int matindex;
	/// marked for refinement
	public int marked;
	/// edge without node k (0,1,2)
	public int markededge;

	public bool incorder;
//C++ TO C# CONVERTER TODO TASK: C# does not allow bit fields:
	public uint order:6;
  }


  public class MarkedIdentification
  {
	// number of points of one face (3 or 4) - or edge (in 2d)
	public int np;
	/// 6 or 8 point numbers - or 4 in 2d
	public PointIndex[] pnums = Arrays.InitializeWithDefaultInstances<PointIndex>(8);
	/// marked for refinement
	public int marked;
	/// edge starting with node k (0,1,2, or 3)
	public int markededge;

	public bool incorder;
//C++ TO C# CONVERTER TODO TASK: C# does not allow bit fields:
	public uint order:6;
  }





  public class MarkedTri
  {
//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//	MarkedTri() = default;
//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//	MarkedTri(const MarkedTri&) = default;
//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//	MarkedTri(MarkedTri &&) = default;
//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//	MarkedTri & operator = (const MarkedTri&) = default;
//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//	MarkedTri & operator = (MarkedTri&&) = default;

	/// three point numbers
	public PointIndex[] pnums = Arrays.InitializeWithDefaultInstances<PointIndex>(3);
	/// three geominfos
	public PointGeomInfo[] pgeominfo = Arrays.InitializeWithDefaultInstances<PointGeomInfo>(3);
	/// marked for refinement
	public int marked;
	/// edge without node k
	public int markededge;
	/// surface id
	public int surfid;

	public bool incorder;
//C++ TO C# CONVERTER TODO TASK: C# does not allow bit fields:
	public uint order:6;
  }



  public class MarkedQuad
  {
	/// point numbers
	public PointIndex[] pnums = Arrays.InitializeWithDefaultInstances<PointIndex>(4);
	///
	public PointGeomInfo[] pgeominfo = Arrays.InitializeWithDefaultInstances<PointGeomInfo>(4);
	/// marked for refinement
	public int marked;
	/// marked edge: 0/2 = vertical, 1/3 = horizontal
	public int markededge;
	/// surface id
	public int surfid;

	public bool incorder;
//C++ TO C# CONVERTER TODO TASK: C# does not allow bit fields:
	public uint order:6;
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void Refinement::Bisect(Mesh & mesh, BisectionOptions & opt, Array<double> * quality_loss) const
}
