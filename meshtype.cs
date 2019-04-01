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


/**************************************************************************/
/* File:   meshtype.hpp                                                   */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Okt. 95                                                    */
/**************************************************************************/

namespace netgen
{

  /*
    Classes for NETGEN
  */



  public enum ELEMENT_TYPE : byte
  {
	SEGMENT = 1,
	SEGMENT3 = 2,
	TRIG = 10,
	QUAD = 11,
	TRIG6 = 12,
	QUAD6 = 13,
	QUAD8 = 14,
	TET = 20,
	TET10 = 21,
	PYRAMID = 22,
	PRISM = 23,
	PRISM12 = 24,
	PRISM15 = 27,
	PYRAMID13 = 28,
	HEX = 25,
	HEX20 = 26
  }

  /*
  typedef int ELEMENT_EDGE[2];      // initial point, end point
  typedef int ELEMENT_FACE[4];      // points, last one is -1 for trig
  */

  public class ELEMENT_EDGE
  {
	public int[] vals = new int[2];
	public int this [uint i]
	{
		get
		{
			return vals[i];
		}
		set
		{
			vals[i] = value;
		}
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int operator [] (uint i) const
	public int this [uint i]
	{
		get
		{
			return vals[i];
		}
	}
  }

  public class ELEMENT_FACE
  {
	public int[] vals = new int[4];
	public int this [uint i]
	{
		get
		{
			return vals[i];
		}
		set
		{
			vals[i] = value;
		}
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int operator [] (uint i) const
	public int this [uint i]
	{
		get
		{
			return vals[i];
		}
	}
  }




  public enum POINTTYPE
  {
	  FIXEDPOINT = 1,
	  EDGEPOINT = 2,
	  SURFACEPOINT = 3,
	  INNERPOINT = 4
  }
  public enum ELEMENTTYPE
  {
	  FREEELEMENT,
	  FIXEDELEMENT
  }
  public enum OPTIMIZEGOAL
  {
	  OPT_QUALITY,
	  OPT_CONFORM,
	  OPT_REST,
	  OPT_WORSTCASE,
	  OPT_LEGAL
  }

  /*
  extern DLL_HEADER int GetTimeStamp();
  extern DLL_HEADER int NextTimeStamp();
  */
  public class PointGeomInfo
  {
	public int trignum; // for STL Meshing
	public double u; // for OCC Meshing
	public double v;

//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//	PointGeomInfo() = default;
	// : trignum(-1), u(0), v(0) { ; }
//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//	PointGeomInfo(const PointGeomInfo&) = default;
//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//	PointGeomInfo(PointGeomInfo &&) = default;
//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//	PointGeomInfo & operator = (const PointGeomInfo&) = default;
//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//	PointGeomInfo & operator = (PointGeomInfo&&) = default;

  }



  public class MultiPointGeomInfo
  {
	private int cnt;
	private PointGeomInfo[] mgi = Arrays.InitializeWithDefaultInstances<PointGeomInfo>(DefineConstants.MULTIPOINTGEOMINFO_MAX);
	public MultiPointGeomInfo()
	{
		cnt = 0;
	}
	public int AddPointGeomInfo(PointGeomInfo gi)
	{
	  for (int k = 0; k < cnt; k++)
	  {
		if (mgi[k].trignum == gi.trignum)
		{
	  return 0;
		}
	  }

	  if (cnt < DefineConstants.MULTIPOINTGEOMINFO_MAX)
	  {
	  mgi[cnt] = gi;
	  cnt++;
	  return 0;
	  }

	  throw new Exception("Please report error: MPGI Size too small\n");
	}

	public void Init()
	{
		cnt = 0;
	}
	public void DeleteAll()
	{
		cnt = 0;
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetNPGI() const
	public int GetNPGI()
	{
		return cnt;
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const PointGeomInfo & GetPGI(int i) const
	public PointGeomInfo GetPGI(int i)
	{
		return new netgen.PointGeomInfo(mgi[i - 1]);
	}
  }


  public class EdgePointGeomInfo
  {
	public int edgenr;
	public int body; // for ACIS
	public double dist; // for 2d meshing
	public double u; // for OCC Meshing
	public double v;

	public EdgePointGeomInfo()
	{
		this.edgenr = 0;
		this.body = 0;
		this.dist = 0.0;
		this.u = 0.0;
		this.v = 0.0;
		;
	}


//C++ TO C# CONVERTER NOTE: This 'CopyFrom' method was converted from the original copy assignment operator:
//ORIGINAL LINE: EdgePointGeomInfo & operator = (const EdgePointGeomInfo & gi2)
	public EdgePointGeomInfo CopyFrom (EdgePointGeomInfo gi2)
	{
	  edgenr = gi2.edgenr;
	  body = gi2.body;
	  dist = gi2.dist;
	  u = gi2.u;
	  v = gi2.v;
	  return this;
	}
  }





  public class PointIndex
  {
	private int i;
//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//	PointIndex() = default;
//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//	PointIndex(const PointIndex&) = default;
//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//	PointIndex(PointIndex &&) = default;
//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//	PointIndex & operator = (const PointIndex&) = default;
//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//	PointIndex & operator = (PointIndex&&) = default;

	public PointIndex(int ai)
	{
		this.i = ai;
		;
	}
	// PointIndex & operator= (const PointIndex &ai) { i = ai.i; return *this; }
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: operator int() const
	public static implicit operator int(PointIndex ImpliedObject)
	{
		return ImpliedObject.i;
	}
	public static PointIndex operator ++(int UnnamedParameter)
	{
		PointIndex hi = new PointIndex(this);
		i++;
		return new netgen.PointIndex(hi);
	}
	public static PointIndex operator --(int UnnamedParameter)
	{
		PointIndex hi = new PointIndex(this);
		i--;
		return new netgen.PointIndex(hi);
	}
	public static PointIndex operator ++(PointIndex ImpliedObject)
	{
		ImpliedObject.i++;
		return ImpliedObject;
	}
	public static PointIndex operator --(PointIndex ImpliedObject)
	{
		ImpliedObject.i--;
		return ImpliedObject;
	}
	public void Invalidate()
	{
		i = PointIndex.BASE-1;
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool IsValid() const
	public bool IsValid()
	{
		return i != PointIndex.BASE-1;
	}
#if BASE0
//C++ TO C# CONVERTER NOTE: Enums must be named in C#, so the following enum has been named by the converter:
	public enum AnonymousEnum
	{
		BASE = 0
	}
#else
//C++ TO C# CONVERTER NOTE: Enums must be named in C#, so the following enum has been named by the converter:
	public enum AnonymousEnum2
	{
		BASE = 1
	}
#endif

	public void DoArchive(Archive ar)
	{
		ar i;
	}
  }
  public class PointIndices < 2> <int N>: INDEX_2
  {
//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//	PointIndices() = default;
	public PointIndices(INDEX_2 i2) : base(i2)
	{
		;
	}
	public PointIndices(PointIndex i1, PointIndex i2) : base(i1, i2)
	{
		;
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: PointIndex operator [] (int i) const
	public new PointIndex this [int i]
	{
		get
		{
			return new PointIndex(base[i]);
		}
	}
	public new PointIndex this [int i]
	{
		get
		{
//C++ TO C# CONVERTER TODO TASK: There is no equivalent to 'reinterpret_cast' in C#:
			return reinterpret_cast<PointIndex&>(base[i]);
		}
		set
		{
//C++ TO C# CONVERTER TODO TASK: There is no equivalent to 'reinterpret_cast' in C#:
			reinterpret_cast<PointIndex&>(base[i]) = value;
		}
	}
	public static PointIndices Sort(PointIndex i1, PointIndex i2)
	{
		return new netgen.INDEX_2(INDEX_2.Sort(new netgen.PointIndex(i1), new netgen.PointIndex(i2)));
	}
  }



  public class ElementIndex
  {
	private int i;
	public ElementIndex()
	{
		;
	}
	public ElementIndex(int ai)
	{
		this.i = ai;
		;
	}
//C++ TO C# CONVERTER NOTE: This 'CopyFrom' method was converted from the original copy assignment operator:
//ORIGINAL LINE: ElementIndex & operator = (const ElementIndex & ai)
	public ElementIndex CopyFrom (ElementIndex ai)
	{
		i = ai.i;
		return this;
	}
//C++ TO C# CONVERTER NOTE: This 'CopyFrom' method was converted from the original copy assignment operator:
//ORIGINAL LINE: ElementIndex & operator = (int ai)
	public ElementIndex CopyFrom (int ai)
	{
		i = ai;
		return this;
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: operator int() const
	public static implicit operator int(ElementIndex ImpliedObject)
	{
		return ImpliedObject.i;
	}
	public static ElementIndex operator ++(int UnnamedParameter)
	{
		return new ElementIndex(i++);
	}
	public static ElementIndex operator --(int UnnamedParameter)
	{
		return new ElementIndex(i--);
	}
	public static ElementIndex operator ++(ElementIndex ImpliedObject)
	{
		++ImpliedObject.i;
		return ImpliedObject;
	}
	public static ElementIndex operator --(ElementIndex ImpliedObject)
	{
		--ImpliedObject.i;
		return ImpliedObject;
	}
  }


  public class SurfaceElementIndex
  {
	private int i;
//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//	SurfaceElementIndex() = default;
	public SurfaceElementIndex(int ai)
	{
		this.i = ai;
		;
	}
	/*
	SurfaceElementIndex & operator= (const SurfaceElementIndex & ai) 
	{ i = ai.i; return *this; }
	*/
//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//	SurfaceElementIndex & operator = (const SurfaceElementIndex & ai) = default;
//C++ TO C# CONVERTER NOTE: This 'CopyFrom' method was converted from the original copy assignment operator:
//ORIGINAL LINE: SurfaceElementIndex & operator = (int ai)
	public SurfaceElementIndex CopyFrom (int ai)
	{
		i = ai;
		return this;
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: operator int() const
	public static implicit operator int(SurfaceElementIndex ImpliedObject)
	{
		return ImpliedObject.i;
	}
	public static SurfaceElementIndex operator ++(int UnnamedParameter)
	{
		SurfaceElementIndex hi = new SurfaceElementIndex(this);
		i++;
//C++ TO C# CONVERTER TODO TASK: The following line was determined to contain a copy constructor call - this should be verified and a copy constructor should be created:
//ORIGINAL LINE: return hi;
		return new netgen.SurfaceElementIndex(hi);
	}
	public static SurfaceElementIndex operator --(int UnnamedParameter)
	{
		SurfaceElementIndex hi = new SurfaceElementIndex(this);
		i--;
//C++ TO C# CONVERTER TODO TASK: The following line was determined to contain a copy constructor call - this should be verified and a copy constructor should be created:
//ORIGINAL LINE: return hi;
		return new netgen.SurfaceElementIndex(hi);
	}
	public static SurfaceElementIndex operator ++(SurfaceElementIndex ImpliedObject)
	{
		++ImpliedObject.i;
		return ImpliedObject;
	}
	public static SurfaceElementIndex operator --(SurfaceElementIndex ImpliedObject)
	{
		--ImpliedObject.i;
		return ImpliedObject;
	}
//C++ TO C# CONVERTER TODO TASK: The += operator cannot be overloaded in C#:
	public static SurfaceElementIndex operator += (int inc)
	{
		i += inc;
		return this;
	}

	public void DoArchive(Archive ar)
	{
		ar i;
	}
  }

  public class SegmentIndex
  {
	private int i;
	public SegmentIndex()
	{
		;
	}
	public SegmentIndex(int ai)
	{
		this.i = ai;
		;
	}
//C++ TO C# CONVERTER NOTE: This 'CopyFrom' method was converted from the original copy assignment operator:
//ORIGINAL LINE: SegmentIndex & operator = (const SegmentIndex & ai)
	public SegmentIndex CopyFrom (SegmentIndex ai)
	{
		i = ai.i;
		return this;
	}
//C++ TO C# CONVERTER NOTE: This 'CopyFrom' method was converted from the original copy assignment operator:
//ORIGINAL LINE: SegmentIndex & operator = (int ai)
	public SegmentIndex CopyFrom (int ai)
	{
		i = ai;
		return this;
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: operator int() const
	public static implicit operator int(SegmentIndex ImpliedObject)
	{
		return ImpliedObject.i;
	}
	public static SegmentIndex operator ++(int UnnamedParameter)
	{
		i++;
		return this;
	}
	public static SegmentIndex operator --(int UnnamedParameter)
	{
		i--;
		return this;
	}
  }




  /**
     Point in the mesh.
     Contains layer (a new feature in 4.3 for overlapping meshes.
  */
  public class MeshPoint : Point < 3>
  {
	private int layer;
	private double singular; // singular factor for hp-refinement
	private POINTTYPE type;


	public MeshPoint()
	{
	  ;
	}

	public MeshPoint(Point < 3> ap, int alayer = 1, POINTTYPE apt = POINTTYPE.INNERPOINT)
	{
		this.Point < 3> = ap.functorMethod;
		this.layer = alayer;
		this.singular = 0.0;
		this.type = new netgen.POINTTYPE(apt);
	  ;
	}

	public void SetPoint(Point < 3> ap)
	{
	  base = ap.functorMethod;
	  layer = 0;
	  singular = 0;
	}

	public void Scale(double factor)
	{
		*testout << "before: " << x[0] << "\n";
		x[0] *= factor;
		x[1] *= factor;
		x[2] *= factor;
		*testout << "after: " << x[0] << "\n";
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetLayer() const
	public int GetLayer()
	{
		return layer;
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: POINTTYPE Type() const
	public POINTTYPE Type()
	{
		return type;
	}
	public void SetType(POINTTYPE at)
	{
		type = at;
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: double Singularity() const
	public double Singularity()
	{
		return singular;
	}
	public void Singularity(double s)
	{
		singular = s;
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool IsSingular() const
	public bool IsSingular()
	{
		return (singular != 0.0);
	}

#if PARALLEL
#if PARALLEL
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	private static int MyGetMPIType_type = 0;
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	private static int MyGetMPIType_htype = 0;

	public static int MyGetMPIType()
	{
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	  static int type = null;
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	  static int htype = null;
	  if (MyGetMPIType_type == 0)
	  {
	  MeshPoint hp = new MeshPoint();
	  int[] blocklen = {3, 1, 1};
	  MPI_Aint[] displ = {(string) hp.x[0] - (string) hp, (string) hp.layer - (string) hp, (string) hp.singular - (string) hp};
	  int[] types = {MPI_DOUBLE, MPI_INT, MPI_DOUBLE};
	  *testout << "displ = " << displ[0] << ", " << displ[1] << ", " << displ[2] << "\n";
	  *testout << "sizeof = " << sizeof(MeshPoint) << "\n";
	  MPI_Type_create_struct(3, blocklen, displ, types, MyGetMPIType_htype);
	  MPI_Type_commit(MyGetMPIType_htype);
	  MPI_Aint lb = new MPI_Aint();
	  MPI_Aint ext = new MPI_Aint();
	  MPI_Type_get_extent(MyGetMPIType_htype, lb, ext);
	  *testout << "lb = " << lb << "\n";
	  *testout << "ext = " << ext << "\n";
	  ext = sizeof(MeshPoint);
	  MPI_Type_create_resized(MyGetMPIType_htype, lb, ext, MyGetMPIType_type);
	  MPI_Type_commit(MyGetMPIType_type);

	  }
	  return (int)MyGetMPIType_type;
	}
#endif

#endif

	public void DoArchive(Archive ar)
	{
	  ar[] x[0] x[1] x & layer & singular = Arrays.InitializeWithDefaultInstances<ar>(2);
	  ar & (byte)(type);
	}
  }







  /**
     Triangle element for surface mesh generation.
  */
  public class Element2d
  {
	/// point numbers
	private PointIndex[] pnum = Arrays.InitializeWithDefaultInstances<PointIndex>(DefineConstants.ELEMENT2D_MAXPOINTS);
	/// geom info of points
	private PointGeomInfo[] geominfo = Arrays.InitializeWithDefaultInstances<PointGeomInfo>(DefineConstants.ELEMENT2D_MAXPOINTS);

	/// surface nr
	private short index;
	///
	private ELEMENT_TYPE typ;
	/// number of points
//C++ TO C# CONVERTER TODO TASK: C# does not allow bit fields:
	private uint np:4;
//C++ TO C# CONVERTER TODO TASK: C# does not allow bit fields:
	private bool badel:1;
//C++ TO C# CONVERTER TODO TASK: C# does not allow bit fields:
	private bool refflag:1; // marked for refinement
//C++ TO C# CONVERTER TODO TASK: C# does not allow bit fields:
	private bool strongrefflag:1;
//C++ TO C# CONVERTER TODO TASK: C# does not allow bit fields:
	private bool deleted:1; // element is deleted

	// Philippose - 08 August 2010
	// Set a new property for each element, to 
	// control whether it is visible or not
//C++ TO C# CONVERTER TODO TASK: C# does not allow bit fields:
	private bool visible:1; // element visible
//C++ TO C# CONVERTER TODO TASK: C# does not allow bit fields:
	private bool is_curved:1; // element is (high order) curved
	/// order for hp-FEM
//C++ TO C# CONVERTER TODO TASK: C# does not allow bit fields:
	private uint orderx:6;
//C++ TO C# CONVERTER TODO TASK: C# does not allow bit fields:
	private uint ordery:6;

	// #ifdef PARALLEL
	// int partitionNumber; 
	// #endif

	/// a linked list for all segments in the same face
	private SurfaceElementIndex next = new SurfaceElementIndex();

	///
//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//	Element2d() = default;
//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//	Element2d(const Element2d &) = default;
//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//	Element2d(Element2d &&) = default;
//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//	Element2d & operator = (const Element2d &) = default;
//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//	Element2d & operator = (Element2d &&) = default;
	///
	/*
	Element2d :: Element2d ()
	{
	  for (int i = 0; i < ELEMENT2D_MAXPOINTS; i++)
	    {
	  pnum[i] = 0;
	  geominfo[i].trignum = 0;
	    }
	  np = 3;
	  index = 0;
	  badel = 0;
	  deleted = 0;
	  visible = 1;
	  typ = TRIG;
	  orderx = ordery = 1;
	  refflag = 1;
	  strongrefflag = false;
	  is_curved = false;
	} 
	*/
	public Element2d(int anp)
	{
	  for (int i = 0; i < DefineConstants.ELEMENT2D_MAXPOINTS; i++)
	  {
	  pnum[i] = 0;
	  geominfo[i].trignum = 0;
	  }
	  np = anp;
	  index = 0;
	  badel = 0;
	  deleted = 0;
	  visible = 1;
	  switch (np)
	  {
		case 3:
			typ = ELEMENT_TYPE.TRIG;
			break;
		case 4:
			typ = ELEMENT_TYPE.QUAD;
			break;
		case 6:
			typ = ELEMENT_TYPE.TRIG6;
			break;
		case 8:
			typ = ELEMENT_TYPE.QUAD8;
			break;
	  }
	  orderx = ordery = 1;
	  refflag = 1;
	  strongrefflag = false;
	  is_curved = (np >= 4); // false;
	}

	///
	public Element2d(ELEMENT_TYPE atyp)
	{
	  for (int i = 0; i < DefineConstants.ELEMENT2D_MAXPOINTS; i++)
	  {
	  pnum[i] = 0;
	  geominfo[i].trignum = 0;
	  }

	  SetType(atyp);

	  index = 0;
	  badel = 0;
	  deleted = 0;
	  visible = 1;
	  orderx = ordery = 1;
	  refflag = 1;
	  strongrefflag = false;
	  is_curved = (np >= 4); // false;
	}

	///
	public Element2d(int pi1, int pi2, int pi3)
	{
	  pnum[0] = pi1;
	  pnum[1] = pi2;
	  pnum[2] = pi3;
	  np = 3;
	  typ = ELEMENT_TYPE.TRIG;
	  pnum[3] = 0;
	  pnum[4] = 0;
	  pnum[5] = 0;

	  for (int i = 0; i < DefineConstants.ELEMENT2D_MAXPOINTS; i++)
	  {
		geominfo[i].trignum = 0;
	  }
	  index = 0;
	  badel = 0;
	  refflag = 1;
	  strongrefflag = false;
	  deleted = 0;
	  visible = 1;
	  orderx = ordery = 1;
	  is_curved = false;
	}

	///
	public Element2d(int pi1, int pi2, int pi3, int pi4)
	{
	  pnum[0] = pi1;
	  pnum[1] = pi2;
	  pnum[2] = pi3;
	  pnum[3] = pi4;
	  np = 4;
	  typ = ELEMENT_TYPE.QUAD;

	  pnum[4] = 0;
	  pnum[5] = 0;

	  for (int i = 0; i < DefineConstants.ELEMENT2D_MAXPOINTS; i++)
	  {
		geominfo[i].trignum = 0;
	  }
	  index = 0;
	  badel = 0;
	  refflag = 1;
	  strongrefflag = false;
	  deleted = 0;
	  visible = 1;
	  orderx = ordery = 1;
	  is_curved = true;
	}

	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: ELEMENT_TYPE GetType() const
	public ELEMENT_TYPE GetType()
	{
		return typ;
	}
	/// 
	public void SetType(ELEMENT_TYPE atyp)
	{
	  typ = atyp;
	  switch (typ)
	  {
	case ELEMENT_TYPE.TRIG:
		np = 3;
		break;
	case ELEMENT_TYPE.QUAD:
		np = 4;
		break;
	case ELEMENT_TYPE.TRIG6:
		np = 6;
		break;
	case ELEMENT_TYPE.QUAD6:
		np = 6;
		break;
	case ELEMENT_TYPE.QUAD8:
		np = 8;
		break;
	default:
	  PrintSysError("Element2d::SetType, illegal type ", (int)typ);
	  break;
	  }
	  is_curved = (np >= 4);
	}
	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetNP() const
	public int GetNP()
	{
		return np;
	}
	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetNV() const
	public int GetNV()
	{
	  if (typ == ELEMENT_TYPE.TRIG || typ == ELEMENT_TYPE.TRIG6)
	  {
		return 3;
	  }
	  else
	  {
#if DEBUG
		  if (typ != ELEMENT_TYPE.QUAD && typ != ELEMENT_TYPE.QUAD6 && typ != ELEMENT_TYPE.QUAD8)
		  {
			PrintSysError("element2d::GetNV not implemented for typ", typ);
		  }
#endif
		  return 4;
	  }
	  /*
	  switch (typ)
	{
	case TRIG:
	case TRIG6: return 3;
	      
	case QUAD:
	case QUAD8:
	case QUAD6: return 4;
	default:
#ifdef DEBUG
	  PrintSysError ("element2d::GetNV not implemented for typ", typ)
#endif
		;
	}
	  return np;
	  */
	}

	///
	public PointIndex this [int i]
	{
		get
		{
			return new netgen.PointIndex(pnum[i]);
		}
		set
		{
			pnum[i] = value;
		}
	}
	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const PointIndex & operator [] (int i) const
	public PointIndex this [int i]
	{
		get
		{
			return new netgen.PointIndex(pnum[i]);
		}
		set
		{
			pnum[i] = value;
		}
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: FlatArray<const PointIndex> PNums() const
	public FlatArray< PointIndex> PNums()
	{
		return new FlatArray<const PointIndex> (np, pnum[0]);
	}
	public FlatArray<PointIndex> PNums()
	{
		return new FlatArray<PointIndex> (np, pnum[0]);
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: auto Vertices() const
//C++ TO C# CONVERTER TODO TASK: The return type of the following function could not be determined:
	public auto Vertices()
	{
		return new FlatArray<const PointIndex> ((uint)GetNV(), pnum[0]);
	}

	///
	public PointIndex PNum(int i)
	{
		return new netgen.PointIndex(pnum[i - 1]);
	}
	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const PointIndex & PNum(int i) const
	public PointIndex PNum(int i)
	{
		return new netgen.PointIndex(pnum[i - 1]);
	}
	///
	public PointIndex PNumMod(int i)
	{
		return new netgen.PointIndex(pnum[(i - 1) % np]);
	}
	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const PointIndex & PNumMod(int i) const
	public PointIndex PNumMod(int i)
	{
		return new netgen.PointIndex(pnum[(i - 1) % np]);
	}
	///

	///
	public PointGeomInfo GeomInfoPi(int i)
	{
		return new netgen.PointGeomInfo(geominfo[i - 1]);
	}
	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const PointGeomInfo & GeomInfoPi(int i) const
	public PointGeomInfo GeomInfoPi(int i)
	{
		return new netgen.PointGeomInfo(geominfo[i - 1]);
	}
	///
	public PointGeomInfo GeomInfoPiMod(int i)
	{
		return new netgen.PointGeomInfo(geominfo[(i - 1) % np]);
	}
	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const PointGeomInfo & GeomInfoPiMod(int i) const
	public PointGeomInfo GeomInfoPiMod(int i)
	{
		return new netgen.PointGeomInfo(geominfo[(i - 1) % np]);
	}

	public void DoArchive(Archive ar)
	{
	  short _np;
	  short _typ;
	  bool _curved;
	  bool _vis;
	  bool _deleted;
	  if (ar.Output())
	  {
			_np = np;
			_typ = (short)typ;
			_curved = is_curved;
		  _vis = visible;
		  _deleted = deleted;
	  }
	  ar & _np & _typ & index & _curved & _vis & _deleted;
	  // ar & next; don't need 
	  if (ar.Input())
	  {
			np = _np;
			typ = ELEMENT_TYPE(_typ);
			is_curved = _curved;
		  visible = _vis;
		  deleted = _deleted;
	  }
	  for (uint i = 0; i < np; i++)
	  {
		ar[] pnum = Arrays.InitializeWithDefaultInstances<ar>(i);
	  }
	}

	public void SetIndex(int si)
	{
		index = (short)si;
	}
	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetIndex() const
	public int GetIndex()
	{
		return index;
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetOrder() const
	public int GetOrder()
	{
		return orderx;
	}
	public void SetOrder(int aorder)
	{
		orderx = ordery = aorder;
	}


//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetOrder(int & ox, int & oy) const
	public void GetOrder(ref int ox, ref int oy)
	{
		ox = orderx, oy = ordery;
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetOrder(int & ox, int & oy, int & oz) const
	public void GetOrder(ref int ox, ref int oy, ref int oz)
	{
		ox = orderx;
		oy = ordery;
		oz = 0;
	}
	public void SetOrder(int ox, int oy, int UnnamedParameter)
	{
		orderx = ox;
		ordery = oy;
	}
	public void SetOrder(int ox, int oy)
	{
		orderx = ox;
		ordery = oy;
	}


	///

	/*
	  void Element2d :: SetType (ELEMENT_TYPE atyp)
	  {
	  typ = atyp;
	  switch (typ)
	  {
	  case TRIG: np = 3; break;
	  case QUAD: np = 4; break;
	  case TRIG6: np = 6; break;
	  case QUAD6: np = 6; break;
	  default:
	  PrintSysError ("Element2d::SetType, illegal type ", typ);
	  }
	  }
	*/


//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetBox(const T_POINTS & points, Box3d & box) const
	public void GetBox(T_POINTS points, Box3d box)
	{
	  box.SetPoint(points.Get(pnum[0]));
	  for (uint i = 1; i < np; i++)
	  {
		box.AddPoint(points.Get(pnum[i]));
	  }
	}

	/// invert orientation
	public void Invert()
	{
	  if (typ == ELEMENT_TYPE.TRIG)
	  {
		netgen.GlobalMembers.Swap(ref PNum(2), ref PNum(3));
	  }
	  else
	  {
		Invert2();
	  }
	}

	///
	public void Invert2()
	{
	  switch (typ)
	  {
		case ELEMENT_TYPE.TRIG:
		{
			netgen.GlobalMembers.Swap(ref pnum[1], ref pnum[2]);
			break;
		}
		case ELEMENT_TYPE.TRIG6:
		{
			netgen.GlobalMembers.Swap(ref pnum[1], ref pnum[2]);
			netgen.GlobalMembers.Swap(ref pnum[4], ref pnum[5]);
			break;
		}
		case ELEMENT_TYPE.QUAD:
		{
			netgen.GlobalMembers.Swap(ref pnum[0], ref pnum[3]);
			netgen.GlobalMembers.Swap(ref pnum[1], ref pnum[2]);
			break;
		}
		default:
		{
			cerr << "Element2d::Invert2, illegal element type " << (int)typ << "\n";
		}
		  break;
	  }
	}

	/// first point number is smallest
	public void NormalizeNumbering()
	{
	  if (GetNP() == 3)
	  {
	  if (PNum(1) < PNum(2) && PNum(1) < PNum(3))
	  {
		return;
	  }
	  else
	  {
		  if (PNum(2) < PNum(3))
		  {
		  PointIndex pi1 = PNum(2);
		  PNum(2) = PNum(3);
		  PNum(3) = PNum(1);
		  PNum(1) = pi1;
		  }
		  else
		  {
		  PointIndex pi1 = PNum(3);
		  PNum(3) = PNum(2);
		  PNum(2) = PNum(1);
		  PNum(1) = pi1;
		  }
	  }
	  }
	  else
	  {
		NormalizeNumbering2();
	  }
	}

	///
	public void NormalizeNumbering2()
	{
	  if (GetNP() == 3)
	  {
		  if (PNum(1) < PNum(2) && PNum(1) < PNum(3))
		  {
			return;
		  }
		  else
		  {
			  if (PNum(2) < PNum(3))
			  {
				  PointIndex pi1 = PNum(2);
				  PNum(2) = PNum(3);
				  PNum(3) = PNum(1);
				  PNum(1) = pi1;
			  }
			  else
			  {
				  PointIndex pi1 = PNum(3);
				  PNum(3) = PNum(2);
				  PNum(2) = PNum(1);
				  PNum(1) = pi1;
			  }
		  }
	  }
	  else
	  {
		  int mini = 1;
		  for (int i = 2; i <= GetNP(); i++)
		  {
			if (PNum(i) < PNum(mini))
			{
				mini = i;
			}
		  }

		  Element2d hel = this;
		  for (int i = 1; i <= GetNP(); i++)
		  {
			PNum(i) = hel.PNumMod(i + mini - 1);
		  }
	  }
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool BadElement() const
	public bool BadElement()
	{
		return badel;
	}

	// friend ostream & operator<<(ostream  & s, const Element2d & el);
//C++ TO C# CONVERTER TODO TASK: C# has no concept of a 'friend' class:
//	friend class Mesh;


	/// get number of 'integration points'
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetNIP() const
	public int GetNIP()
	{
	  int nip;
	  switch (np)
	  {
		case 3:
			nip = 1;
			break;
		case 4:
			nip = 4;
			break;
		default:
			nip = 0;
			break;
	  }
	  return nip;
	}

//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	private double[][] GetIntegrationPoint_eltriqp =
	{
		new double[] {1.0 / 3.0, 1.0 / 3.0, 0.5}
	};
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	private double[][] GetIntegrationPoint_elquadqp =
	{
		new double[] {0, 0, 0.25},
		new double[] {0, 1, 0.25},
		new double[] {1, 0, 0.25},
		new double[] {1, 1, 0.25}
	};

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetIntegrationPoint(int ip, Point2d & p, double & weight) const
	public void GetIntegrationPoint(int ip, Point2d p, ref double weight)
	{
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	  static double eltriqp[1][3] = { { 1.0/3.0, 1.0/3.0, 0.5 } };

//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	  static double elquadqp[4][3] = { { 0, 0, 0.25 }, { 0, 1, 0.25 }, { 1, 0, 0.25 }, { 1, 1, 0.25 } };

//C++ TO C# CONVERTER TODO TASK: C# does not have an equivalent to pointers to value types:
//ORIGINAL LINE: double * pp = 0;
	  double pp = null;
	  switch (typ)
	  {
		case ELEMENT_TYPE.TRIG:
			pp = GetIntegrationPoint_eltriqp[0][0];
			break;
		case ELEMENT_TYPE.QUAD:
			pp = GetIntegrationPoint_elquadqp[ip - 1][0];
			break;
		default:
		  PrintSysError("Element2d::GetIntegrationPoint, illegal type ", (int)typ);
		  break;
	  }

	  p.X() = pp[0];
	  p.Y() = pp[1];
	  weight = pp[2];
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetTransformation(int ip, const Array<Point2d> & points, class DenseMatrix & trans) const;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//	void GetTransformation(int ip, Array<Point2d> points, DenseMatrix trans);
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetTransformation(int ip, class DenseMatrix & pmat, class DenseMatrix & trans) const
	public void GetTransformation(int ip, DenseMatrix pmat, DenseMatrix trans)
	{
	  //  int np = GetNP();

#if DEBUG
	  if (pmat.Width() != np || pmat.Height() != 2)
	  {
		  (*testout) << "GetTransofrmation: pmat doesn't fit" << "\n";
		  return;
	  }
#endif

	  ComputeIntegrationPointData();
	  DenseMatrix dshapep = null;
	  switch (typ)
	  {
		case ELEMENT_TYPE.TRIG:
			dshapep = GlobalMembers.ipdtrig.Get(ip).dshape;
			break;
		case ELEMENT_TYPE.QUAD:
			dshapep = GlobalMembers.ipdquad.Get(ip).dshape;
			break;
		default:
		  PrintSysError("Element2d::GetTransformation, illegal type ", (int)typ);
		  break;
	  }

	  CalcABt(pmat.functorMethod, dshapep.functorMethod, trans.functorMethod);
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetShape(const Point2d & p, class Vector & shape) const;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//	void GetShape(Point2d p, Vector shape);
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetShapeNew(const Point<2> & p, class FlatVector & shape) const;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//	void GetShapeNew(Point<2> p, FlatVector shape);
//C++ TO C# CONVERTER TODO TASK: The original C++ template specifier was replaced with a C# generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <typename T>
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetShapeNew(const Point<2,T> & p, TFlatVector<T> shape) const
	public void GetShapeNew<T>(Point<2,T> p, TFlatVector<T> shape)
	{
	  switch (typ)
	  {
		case ELEMENT_TYPE.TRIG:
		{
			shape.functorMethod(0) = p(0);
			shape.functorMethod(1) = p(1);
			shape.functorMethod(2) = 1 - p(0) - p(1);
			break;
		}

		case ELEMENT_TYPE.QUAD:
		{
			shape.functorMethod(0) = (1 - p(0)) * (1 - p(1));
			shape.functorMethod(1) = p(0) * (1 - p(1));
			shape.functorMethod(2) = p(0) * p(1);
			shape.functorMethod(3) = (1 - p(0)) * p(1);
			break;
		}
		default:
		  throw new Exception("illegal element type in GetShapeNew");
	  }
	}

	/// matrix 2 * np
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetDShape(const Point2d & p, class DenseMatrix & dshape) const;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//	void GetDShape(Point2d p, DenseMatrix dshape);
//C++ TO C# CONVERTER TODO TASK: The original C++ template specifier was replaced with a C# generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <typename T>
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetDShapeNew(const Point<2,T> & p, class MatrixFixWidth<2,T> & dshape) const;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//	void GetDShapeNew<T>(Point<2,T> p, MatrixFixWidth<2,T> dshape);

	/// matrix 2 * np
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetPointMatrix(const Array<Point2d> & points, class DenseMatrix & pmat) const;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//	void GetPointMatrix<T>(Array<Point2d> points, DenseMatrix pmat);

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void ComputeIntegrationPointData() const
	public void ComputeIntegrationPointData<T>()
	{
	  switch (np)
	  {
		case 3:
			if (GlobalMembers.ipdtrig.Size())
			{
				return;
			}
			break;
		case 4:
			if (GlobalMembers.ipdquad.Size())
			{
				return;
			}
			break;
	  }

	  for (int i = 1; i <= GetNIP(); i++)
	  {
		  IntegrationPointData ipd = new IntegrationPointData();
		  Point2d hp = new Point2d();
		  GetIntegrationPoint(i, hp, ref ipd.weight);
		  ipd.p(0) = hp.X();
		  ipd.p(1) = hp.Y();
		  ipd.p(2) = 0;

		  ipd.shape.SetSize(GetNP());
		  ipd.dshape.SetSize(2, GetNP());

		  GetShape(hp, ipd.shape);
		  GetDShape(hp, ipd.dshape.functorMethod);

		  switch (np)
		  {
			case 3:
				GlobalMembers.ipdtrig.Append(ipd);
				break;
			case 4:
				GlobalMembers.ipdquad.Append(ipd);
				break;
		  }
	  }
	}


//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: double CalcJacobianBadness(const Array<Point2d> & points) const
	public double CalcJacobianBadness(Array<Point2d> points)
	{
	  int i;
	  int j;
	  int nip = GetNIP();
	  DenseMatrix trans = new DenseMatrix(2, 2);
	  DenseMatrix pmat = new DenseMatrix();

	  pmat.SetSize(2, GetNP());
	  GetPointMatrix(points, pmat.functorMethod);

	  double err = 0;
	  for (i = 1; i <= nip; i++)
	  {
		  GetTransformation(i, pmat.functorMethod, trans.functorMethod);

		  // Frobenius norm
		  double frob = 0;
		  for (j = 1; j <= 4; j++)
		  {
			frob += netgen.GlobalMembers.sqr(trans.Get(j));
		  }
		  frob = ngsimd.GlobalMembers.sqrt(frob);
		  frob /= 2;

		  double det = trans.Det();

		  if (det <= 0)
		  {
			err += 1e12;
		  }
		  else
		  {
			err += frob * frob / det;
		  }
	  }

	  err /= nip;
	  return err;
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: double CalcJacobianBadness(const T_POINTS & points, const Vec<3> & n) const
	public double CalcJacobianBadness(T_POINTS points, Vec < 3> n)
	{
	  int i;
	  int j;
	  int nip = GetNIP();
	  DenseMatrix trans = new DenseMatrix(2, 2);
	  DenseMatrix pmat = new DenseMatrix();

	  pmat.SetSize(2, GetNP());

	  Vec < 3> t1, t2;
	  t1 = n.GetNormal();
	  t2 = netgen.GlobalMembers.Cross(n.functorMethod, t1);

	  for (i = 1; i <= GetNP(); i++)
	  {
		  Point3d p = points.Get(PNum(i));
		  pmat.Elem(1, i) = p.X() * t1(0) + p.Y() * t1(1) + p.Z() * t1(2);
		  pmat.Elem(2, i) = p.X() * t2(0) + p.Y() * t2(1) + p.Z() * t2(2);
	  }

	  double err = 0;
	  for (i = 1; i <= nip; i++)
	  {
		  GetTransformation(i, pmat.functorMethod, trans.functorMethod);

		  // Frobenius norm
		  double frob = 0;
		  for (j = 1; j <= 4; j++)
		  {
			frob += netgen.GlobalMembers.sqr(trans.Get(j));
		  }
		  frob = ngsimd.GlobalMembers.sqrt(frob);
		  frob /= 2;

		  double det = trans.Det();
		  if (det <= 0)
		  {
			err += 1e12;
		  }
		  else
		  {
			err += frob * frob / det;
		  }
	  }

	  err /= nip;
	  return err;
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: double CalcJacobianBadnessDirDeriv(const Array<Point2d> & points, int pi, Vec2d & dir, double & dd) const
	public double CalcJacobianBadnessDirDeriv(Array<Point2d> points, int pi, Vec2d dir, ref double dd)
	{
	  if (typ == ELEMENT_TYPE.QUAD)
	  {
		  Mat < 2,2> trans, dtrans;
		  Mat < 2,4> vmat, pmat;

		  for (int j = 0; j < 4; j++)
		  {
			  Point2d p = points.Get(this[j]);
			  pmat(0, j) = p.X();
			  pmat(1, j) = p.Y();
		  }

		  vmat = 0.0;
		  vmat(0, pi - 1) = dir.X();
		  vmat(1, pi - 1) = dir.Y();

		  double err = 0;
		  dd = 0;

		  for (int i = 0; i < 4; i++)
		  {
			  int ix1 = GlobalMembers.qip_table[i][0];
			  int ix2 = GlobalMembers.qip_table[i][1];
			  int iy1 = GlobalMembers.qip_table[i][2];
			  int iy2 = GlobalMembers.qip_table[i][3];

			  trans(0,0) = pmat(0, ix2) - pmat(0,ix1);
			  trans(1,0) = pmat(1, ix2) - pmat(1,ix1);
			  trans(0,1) = pmat(0, iy2) - pmat(0,iy1);
			  trans(1,1) = pmat(1, iy2) - pmat(1,iy1);

			  double det = trans(0,0) * trans(1,1) - trans(1,0) * trans(0,1);

			  if (det <= 0)
			  {
				  dd = 0;
				  return 1e12;
			  }

			  dtrans(0,0) = vmat(0, ix2) - vmat(0,ix1);
			  dtrans(1,0) = vmat(1, ix2) - vmat(1,ix1);
			  dtrans(0,1) = vmat(0, iy2) - vmat(0,iy1);
			  dtrans(1,1) = vmat(1, iy2) - vmat(1,iy1);


			  // Frobenius norm
			  double frob = 0;
			  for (int j = 0; j < 4; j++)
			  {
				frob += netgen.GlobalMembers.sqr(trans(j));
			  }
			  frob = ngsimd.GlobalMembers.sqrt(frob);

			  double dfrob = 0;
			  for (int j = 0; j < 4; j++)
			  {
				dfrob += trans(j) * dtrans(j);
			  }
			  dfrob = dfrob / frob;

			  frob /= 2;
			  dfrob /= 2;


			  // ddet = \sum_j det (m_j)   with m_j = trans, except col j = dtrans
			  double ddet = dtrans(0,0) * trans(1,1) - trans(0,1) * dtrans(1,0) + trans(0,0) * dtrans(1,1) - dtrans(0,1) * trans(1,0);

			  err += frob * frob / det;
			  dd += (2 * frob * dfrob * det - frob * frob * ddet) / (det * det);
		  }

		  err /= 4;
		  dd /= 4;
		  return err;
	  }

	  int nip = GetNIP();
	  DenseMatrix trans = new DenseMatrix(2, 2);
	  DenseMatrix dtrans = new DenseMatrix(2, 2);
	  DenseMatrix pmat = new DenseMatrix();
	  DenseMatrix vmat = new DenseMatrix();

	  pmat.SetSize(2, GetNP());
	  vmat.SetSize(2, GetNP());

	  GetPointMatrix(points, pmat.functorMethod);

	  vmat = 0.0;
	  vmat.Elem(1, pi) = dir.X();
	  vmat.Elem(2, pi) = dir.Y();


	  double err = 0;
	  dd = 0;

	  for (int i = 1; i <= nip; i++)
	  {
		  GetTransformation(i, pmat.functorMethod, trans.functorMethod);
		  GetTransformation(i, vmat.functorMethod, dtrans.functorMethod);

		  // Frobenius norm
		  double frob = 0;
		  for (int j = 1; j <= 4; j++)
		  {
			frob += netgen.GlobalMembers.sqr(trans.Get(j));
		  }
		  frob = ngsimd.GlobalMembers.sqrt(frob);

		  double dfrob = 0;
		  for (int j = 1; j <= 4; j++)
		  {
			dfrob += trans.Get(j) * dtrans.Get(j);
		  }
		  dfrob = dfrob / frob;

		  frob /= 2;
		  dfrob /= 2;

		  double det = trans.functorMethod(0,0) * trans.functorMethod(1,1) - trans.functorMethod(1,0) * trans.functorMethod(0,1);

		  // ddet = \sum_j det (m_j)   with m_j = trans, except col j = dtrans
		  double ddet = dtrans.functorMethod(0,0) * trans.functorMethod(1,1) - trans.functorMethod(0,1) * dtrans.functorMethod(1,0) + trans.functorMethod(0,0) * dtrans.functorMethod(1,1) - dtrans.functorMethod(0,1) * trans.functorMethod(1,0);

		  if (det <= 0)
		  {
			err += 1e12;
		  }
		  else
		  {
			  err += frob * frob / det;
			  dd += (2 * frob * dfrob * det - frob * frob * ddet) / (det * det);
		  }
	  }

	  err /= nip;
	  dd /= nip;
	  return err;
	}



	public void Delete()
	{
	  deleted = 1;
	  foreach (PointIndex p in pnum)
	  {
		  p.Invalidate();
	  }
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool IsDeleted() const
	public bool IsDeleted()
	{
#if DEBUG
	  if (pnum[0] < PointIndex.BASE && !deleted)
	  {
	cerr << "Surfelement has illegal pnum, but not marked as deleted" << "\n";
	  }
#endif
	  return deleted;
	}

	// Philippose - 08 August 2010
	// Access functions for the new property: visible
	public void Visible(bool vis = true)
	{
		visible = vis;
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool IsVisible() const
	public bool IsVisible()
	{
		return visible;
	}

	public void SetRefinementFlag(bool rflag = true)
	{
		refflag = rflag;
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool TestRefinementFlag() const
	public bool TestRefinementFlag()
	{
		return refflag;
	}

	public void SetStrongRefinementFlag(bool rflag = true)
	{
		strongrefflag = rflag;
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool TestStrongRefinementFlag() const
	public bool TestStrongRefinementFlag()
	{
		return strongrefflag;
	}


//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool IsCurved() const
	public bool IsCurved()
	{
		return is_curved;
	}
	public void SetCurved(bool acurved)
	{
		is_curved = acurved;
	}

	public SurfaceElementIndex NextElement()
	{
//C++ TO C# CONVERTER TODO TASK: The following line was determined to contain a copy constructor call - this should be verified and a copy constructor should be created:
//ORIGINAL LINE: return next;
		return new netgen.SurfaceElementIndex(next);
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool operator ==(const Element2d & el2) const
	public static bool operator == (Element2d ImpliedObject, Element2d el2)
	{
	  bool retval = (el2.GetNP() == np);
	  for (int i = 0; retval && i < np; i++)
	  {
		retval = (el2[i] == (*ImpliedObject)[i]);
	  }

	  return retval;
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int HasFace(const Element2d & el) const
	public int HasFace(Element2d el)
	{
	  //nur für tets!!! hannes
	  for (int i = 1; i <= 3; i++)
	  {
		  if (PNumMod(i) == el[0] && PNumMod(i + 1) == el[1] && PNumMod(i + 2) == el[2])
		  {
			  return 1;
		  }
	  }
	  return 0;
	}

	///
	public int meshdocval;
	///
	public int hp_elnr;

	/*
#ifdef PARALLEL
	int GetPartition () const { return partitionNumber; }
	void SetPartition (int nr) { partitionNumber = nr; }; 
#endif
	*/
  }





  public class IntegrationPointData
  {
	public Point < 3> p;
	public double weight;
	public Vector shape = new Vector();
	public DenseMatrix dshape = new DenseMatrix();
  }








  /**
     Volume element
  */
  public class Element
  {
	/// point numbers
	private PointIndex[] pnum = Arrays.InitializeWithDefaultInstances<PointIndex>(DefineConstants.ELEMENT_MAXPOINTS);
	///
	private ELEMENT_TYPE typ;
	/// number of points (4..tet, 5..pyramid, 6..prism, 8..hex, 10..quad tet, 12..quad prism)
//C++ TO C# CONVERTER TODO TASK: C# does not allow bit fields:
	private int np:6;
	///
	private class flagstruct
	{
//C++ TO C# CONVERTER TODO TASK: C# does not allow bit fields:
	  public bool marked:1; // marked for refinement
//C++ TO C# CONVERTER TODO TASK: C# does not allow bit fields:
	  public bool badel:1; // angles worse then limit
//C++ TO C# CONVERTER TODO TASK: C# does not allow bit fields:
	  public bool reverse:1; // for refinement a la Bey
//C++ TO C# CONVERTER TODO TASK: C# does not allow bit fields:
	  public bool illegal:1; // illegal, will be split or swapped
//C++ TO C# CONVERTER TODO TASK: C# does not allow bit fields:
	  public bool illegal_valid:1; // is illegal-flag valid ?
//C++ TO C# CONVERTER TODO TASK: C# does not allow bit fields:
	  public bool badness_valid:1; // is badness valid ?
//C++ TO C# CONVERTER TODO TASK: C# does not allow bit fields:
	  public bool refflag:1; // mark element for refinement
//C++ TO C# CONVERTER TODO TASK: C# does not allow bit fields:
	  public bool strongrefflag:1;
//C++ TO C# CONVERTER TODO TASK: C# does not allow bit fields:
	  public bool deleted:1; // element is deleted, will be removed from array
//C++ TO C# CONVERTER TODO TASK: C# does not allow bit fields:
	  public bool @fixed:1; // don't change element in optimization
	}

	/// sub-domain index
	private short index;
	/// order for hp-FEM
//C++ TO C# CONVERTER TODO TASK: C# does not allow bit fields:
	private uint orderx:6;
//C++ TO C# CONVERTER TODO TASK: C# does not allow bit fields:
	private uint ordery:6;
//C++ TO C# CONVERTER TODO TASK: C# does not allow bit fields:
	private uint orderz:6;
	/* unsigned int levelx:6;
	   unsigned int levely:6;
	   unsigned int levelz:6; */ 
	/// stored shape-badness of element
	private float badness;
//C++ TO C# CONVERTER TODO TASK: C# does not allow bit fields:
	private bool is_curved:1; // element is (high order) curved

	// #ifdef PARALLEL
	/// number of partition for parallel computation 
	// int partitionNumber;

	// #endif

	public flagstruct flags = new flagstruct();

	///
//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//	DLL_HEADER Element() = default;
//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//	Element(const Element &) = default;
//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//	Element(Element &&) = default;
//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//	Element & operator = (const Element &) = default;
//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//	Element & operator = (Element &&) = default;

	///

	/*
	Element :: Element ()
	{
	  typ = TET;
	  np = 4;
	  for (int i = 0; i < ELEMENT_MAXPOINTS; i++)
	    pnum[i] = 0;
	  index = 0;
	  flags.marked = 1;
	  flags.badel = 0;
	  flags.reverse = 0;
	  flags.illegal = 0;
	  flags.illegal_valid = 0;
	  flags.badness_valid = 0;
	  flags.refflag = 1;
	  flags.strongrefflag = false;
	  flags.deleted = 0;
	  flags.fixed = 0;
	  orderx = ordery = orderz = 1;
	  is_curved = false;
  #ifdef PARALLEL
	  partitionNumber = -1;
  #endif
	}
	*/

	public Element(int anp)
	{
	  np = anp;
	  int i;
	  for (i = 0; i < DefineConstants.ELEMENT_MAXPOINTS; i++)
	  {
		pnum[i] = 0;
	  }
	  index = 0;
	  flags.marked = 1;
	  flags.badel = 0;
	  flags.reverse = 0;
	  flags.illegal = 0;
	  flags.illegal_valid = 0;
	  flags.badness_valid = 0;
	  flags.refflag = 1;
	  flags.strongrefflag = false;
	  flags.deleted = 0;
	  flags.@fixed = 0;

	  switch (np)
	  {
		case 4:
			typ = ELEMENT_TYPE.TET;
			break;
		case 5:
			typ = ELEMENT_TYPE.PYRAMID;
			break;
		case 6:
			typ = ELEMENT_TYPE.PRISM;
			break;
		case 8:
			typ = ELEMENT_TYPE.HEX;
			break;
		case 10:
			typ = ELEMENT_TYPE.TET10;
			break;
		case 13:
			typ = ELEMENT_TYPE.PYRAMID13;
			break;
		case 15:
			typ = ELEMENT_TYPE.PRISM15;
			break;
		case 20:
			typ = ELEMENT_TYPE.HEX20;
			break;
		default:
			cerr << "Element::Element: unknown element with " << np << " points" << "\n";
			break;
	  }
	  orderx = ordery = orderz = 1;
	  is_curved = typ != ELEMENT_TYPE.TET; // false;
	}

	///
	public Element(ELEMENT_TYPE type)
	{
	  SetType(type);

	  int i;
	  for (i = 0; i < DefineConstants.ELEMENT_MAXPOINTS; i++)
	  {
		pnum[i] = 0;
	  }
	  index = 0;
	  flags.marked = 1;
	  flags.badel = 0;
	  flags.reverse = 0;
	  flags.illegal = 0;
	  flags.illegal_valid = 0;
	  flags.badness_valid = 0;
	  flags.refflag = 1;
	  flags.strongrefflag = false;
	  flags.deleted = 0;
	  flags.@fixed = 0;
	  orderx = ordery = orderz = 1;
	  is_curved = typ != ELEMENT_TYPE.TET; // false;
	  // #ifdef PARALLEL
	  // partitionNumber = -1;
	  // #endif
	}

	///
	// Element & operator= (const Element & el2);

	///

	/*
	Element & Element :: operator= (const Element & el2)
	{
	  typ = el2.typ;
	  np = el2.np;
	  for (int i = 0; i < ELEMENT_MAXPOINTS; i++)
	    pnum[i] = el2.pnum[i];
	  index = el2.index;
	  flags = el2.flags;
	  orderx = el2.orderx;
	  ordery = el2.ordery;
	  orderz = el2.orderz;
	  hp_elnr = el2.hp_elnr;
	  flags = el2.flags;
	  is_curved = el2.is_curved;
	  return *this;
	}
	*/


	public void SetNP(int anp)
	{
	  np = anp;
	  switch (np)
	  {
		case 4:
			typ = ELEMENT_TYPE.TET;
			break;
		case 5:
			typ = ELEMENT_TYPE.PYRAMID;
			break;
		case 6:
			typ = ELEMENT_TYPE.PRISM;
			break;
		case 8:
			typ = ELEMENT_TYPE.HEX;
			break;
		case 10:
			typ = ELEMENT_TYPE.TET10;
			break;
		case 13:
			typ = ELEMENT_TYPE.PYRAMID13;
			break;
		case 15:
			typ = ELEMENT_TYPE.PRISM15;
			break;
		case 20:
			typ = ELEMENT_TYPE.HEX20;
			break;
		  //
		default:
			break;
		  cerr << "Element::SetNP unknown element with " << np << " points" << "\n";
		  break;
	  }
	}

	///
	public void SetType(ELEMENT_TYPE atyp)
	{
	  typ = atyp;
	  switch (atyp)
	  {
		case ELEMENT_TYPE.TET:
			np = 4;
			break;
		case ELEMENT_TYPE.PYRAMID:
			np = 5;
			break;
		case ELEMENT_TYPE.PRISM:
			np = 6;
			break;
		case ELEMENT_TYPE.HEX:
			np = 8;
			break;
		case ELEMENT_TYPE.TET10:
			np = 10;
			break;
		case ELEMENT_TYPE.PYRAMID13:
			np = 13;
			break;
		case ELEMENT_TYPE.PRISM12:
			np = 12;
			break;
		case ELEMENT_TYPE.PRISM15:
			np = 15;
			break;
		case ELEMENT_TYPE.HEX20:
			np = 20;
			break;

		default:
			break;
		  cerr << "Element::SetType unknown type  " << (int)typ << "\n";
		  break;
	  }
	}

	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetNP() const
	public int GetNP()
	{
		return np;
	}
	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: byte GetNV() const
	public byte GetNV()
	{
	  __assume(typ >= ELEMENT_TYPE.TET && typ <= ELEMENT_TYPE.PYRAMID13);
	  switch (typ)
	  {
		case ELEMENT_TYPE.TET:
		case ELEMENT_TYPE.TET10:
		  return 4;
		case ELEMENT_TYPE.PRISM12:
		case ELEMENT_TYPE.PRISM15:
		case ELEMENT_TYPE.PRISM:
	  return 6;
	case ELEMENT_TYPE.PYRAMID:
		case ELEMENT_TYPE.PYRAMID13:
	  return 5;
	case ELEMENT_TYPE.HEX:
	case ELEMENT_TYPE.HEX20:
	  return 8;
		default: // not a 3D element
#if DEBUG
		  PrintSysError("Element3d::GetNV not implemented for typ ", typ);
#endif
		  __assume(false);
		  return -1;
	  }
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool operator ==(const Element & el2) const
	public static bool operator == (Element ImpliedObject, Element el2)
	{
	  bool retval = (el2.GetNP() == np);
	  for (int i = 0; retval && i < np; i++)
	  {
		retval = (el2[i] == (*ImpliedObject)[i]);
	  }

	  return retval;
	}

	// old style:
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int NP() const
	public int NP()
	{
		return np;
	}

	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: ELEMENT_TYPE GetType() const
	public ELEMENT_TYPE GetType()
	{
		return typ;
	}

	///
	public PointIndex this [int i]
	{
		get
		{
			return new netgen.PointIndex(pnum[i]);
		}
		set
		{
			pnum[i] = value;
		}
	}
	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const PointIndex & operator [] (int i) const
	public PointIndex this [int i]
	{
		get
		{
			return new netgen.PointIndex(pnum[i]);
		}
		set
		{
			pnum[i] = value;
		}
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: FlatArray<const PointIndex> PNums() const
	public FlatArray< PointIndex> PNums()
	{
		return new FlatArray<const PointIndex> (np, pnum[0]);
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: FlatArray<const PointIndex> Vertices() const
	public FlatArray< PointIndex> Vertices()
	{
//C++ TO C# CONVERTER TODO TASK: The following line could not be converted:
		return {GetNV(), &pnum[0]};
	}

	///
	public PointIndex PNum(int i)
	{
		return new netgen.PointIndex(pnum[i - 1]);
	}
	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const PointIndex & PNum(int i) const
	public PointIndex PNum(int i)
	{
		return new netgen.PointIndex(pnum[i - 1]);
	}
	///
	public PointIndex PNumMod(int i)
	{
		return new netgen.PointIndex(pnum[(i - 1) % np]);
	}
	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const PointIndex & PNumMod(int i) const
	public PointIndex PNumMod(int i)
	{
		return new netgen.PointIndex(pnum[(i - 1) % np]);
	}

	public void DoArchive(Archive ar)
	{
	  short _np;
	  short _typ;
	  if (ar.Output())
	  {
			_np = np;
			_typ = (short)typ;
	  }
	  ar & _np & _typ & index;
	  if (ar.Input())
	  {
			np = _np;
			typ = ELEMENT_TYPE(_typ);
	  }
	  for (uint i = 0; i < np; i++)
	  {
		ar[] pnum = Arrays.InitializeWithDefaultInstances<ar>(i);
	  }
	}

	///
	public void SetIndex(int si)
	{
		index = (short)si;
	}
	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetIndex() const
	public int GetIndex()
	{
		return index;
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetOrder() const
	public int GetOrder()
	{
		return orderx;
	}
	public void SetOrder(int aorder)
	{
	  orderx = aorder;
	  ordery = aorder;
	  orderz = aorder;
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetOrder(int & ox, int & oy, int & oz) const
	public void GetOrder(ref int ox, ref int oy, ref int oz)
	{
		ox = orderx;
		oy = ordery;
		oz = orderz;
	}
	public void SetOrder(int ox, int oy, int oz)
	{
	  orderx = ox;
	  ordery = oy;
	  orderz = oz;
	}

	// void GetLevel (int & ox, int & oy, int & oz) const { ox = levelx; oy = levely; oz = levelz; }
	// void SetLevel (int ox, int oy, int oz) { levelx = ox; levely = oy; levelz = oz; }


	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetBox(const T_POINTS & points, Box3d & box) const
	public void GetBox(T_POINTS points, Box3d box)
	{
	  box.SetPoint(points.Get(PNum(1)));
	  box.AddPoint(points.Get(PNum(2)));
	  box.AddPoint(points.Get(PNum(3)));
	  box.AddPoint(points.Get(PNum(4)));
	}

	/// Calculates Volume of element
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: double Volume(const T_POINTS & points) const
	public double Volume(T_POINTS points)
	{
	  Vec < 3> v1 = points.Get(PNum(2)) - points.Get(PNum(1));
	  Vec < 3> v2 = points.Get(PNum(3)) - points.Get(PNum(1));
	  Vec < 3> v3 = points.Get(PNum(4)) - points.Get(PNum(1));

	  return -(netgen.GlobalMembers.Cross(v1, v2) * v3) / 6;
	}

	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void Print(ostream & ost) const
	public void Print(ostream ost)
	{
	  ost << np << " Points: ";
	  for (int i = 1; i <= np; i++)
	  {
		ost << pnum[i - 1] << " " << "\n";
	  }
	}

	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetNFaces() const
	public int GetNFaces()
	{
	  switch (typ)
	  {
	case ELEMENT_TYPE.TET:
	case ELEMENT_TYPE.TET10:
		return 4;
	case ELEMENT_TYPE.PYRAMID:
case ELEMENT_TYPE.PYRAMID13:
	return 5;
	case ELEMENT_TYPE.PRISM:
		case ELEMENT_TYPE.PRISM15:
	case ELEMENT_TYPE.PRISM12:
		return 5;
		case ELEMENT_TYPE.HEX:
	case ELEMENT_TYPE.HEX20:
		  return 6;
	default:
#if DEBUG
//C++ TO C# CONVERTER TODO TASK: Statements that are interrupted by preprocessor statements are not converted by C++ to C# Converter:
	  PrintSysError("element3d::GetNFaces not implemented for typ", typ)
#endif
		;
		break;
	  }
	  return 0;
	}
	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: inline void GetFace(int i, Element2d & face) const
	public void GetFace(int i, Element2d face)
	{
	  if (typ == ELEMENT_TYPE.TET)
	  {
	  face.SetType(ELEMENT_TYPE.TRIG);
	  face[0] = pnum[GlobalMembers.gftetfacesa[i - 1][0]];
	  face[1] = pnum[GlobalMembers.gftetfacesa[i - 1][1]];
	  face[2] = pnum[GlobalMembers.gftetfacesa[i - 1][2]];
	  }
	  else
	  {
		GetFace2(i, face);
	  }
	}

	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetFace2(int i, Element2d & face) const
	public void GetFace2(int i, Element2d face)
	{
	  int[][] tetfaces =
	  {
		  new int[] {3, 2, 3, 4, 0},
		  new int[] {3, 3, 1, 4, 0},
		  new int[] {3, 1, 2, 4, 0},
		  new int[] {3, 2, 1, 3, 0}
	  };

	  int[][] tet10faces =
	  {
		  new int[] {3, 2, 3, 4, 10, 9, 8},
		  new int[] {3, 3, 1, 4, 7, 10, 6},
		  new int[] {3, 1, 2, 4, 9, 7, 5},
		  new int[] {3, 2, 1, 3, 6, 8, 5}
	  };

	  int[][] pyramidfaces =
	  {
		  new int[] {4, 1, 4, 3, 2},
		  new int[] {3, 1, 2, 5, 0},
		  new int[] {3, 2, 3, 5, 0},
		  new int[] {3, 3, 4, 5, 0},
		  new int[] {3, 4, 1, 5, 0}
	  };

	  int[][] prismfaces =
	  {
		  new int[] {3, 1, 3, 2, 0},
		  new int[] {3, 4, 5, 6, 0},
		  new int[] {4, 1, 2, 5, 4},
		  new int[] {4, 2, 3, 6, 5},
		  new int[] {4, 3, 1, 4, 6}
	  };

	  int[][] hexfaces =
	  {
		  new int[] {4, 4, 3, 2, 1},
		  new int[] {4, 3, 7, 6, 2},
		  new int[] {4, 7, 8, 5, 6},
		  new int[] {4, 8, 4, 1, 5},
		  new int[] {4, 1, 2, 6, 5},
		  new int[] {4, 3, 4, 8, 7}
	  };


	  switch (np)
	  {
		case 4: // tet
		{
			face.SetType(ELEMENT_TYPE.TRIG);
			for (int j = 1; j <= 3; j++)
			{
			  face.PNum(j) = PNum(tetfaces[i - 1][j]);
			}
			break;
		}

		case 10: // tet10
		{
			face.SetType(ELEMENT_TYPE.TRIG6);
			for (int j = 1; j <= 6; j++)
			{
			  face.PNum(j) = PNum(tet10faces[i - 1][j]);
			}
			break;
		}

		case 5: // pyramid
		{
			// face.SetNP(pyramidfaces[i-1][0]);
			face.SetType((i == 1) ? ELEMENT_TYPE.QUAD : ELEMENT_TYPE.TRIG);
			for (int j = 1; j <= face.GetNP(); j++)
			{
			  face.PNum(j) = PNum(pyramidfaces[i - 1][j]);
			}
			break;
		}
		case 6: // prism
		{
			//	face.SetNP(prismfaces[i-1][0]);
			face.SetType((i >= 3) ? ELEMENT_TYPE.QUAD : ELEMENT_TYPE.TRIG);
			for (int j = 1; j <= face.GetNP(); j++)
			{
			  face.PNum(j) = PNum(prismfaces[i - 1][j]);
			}
			break;
		}
		case 8:
		{
			face.SetType(ELEMENT_TYPE.QUAD);
			for (int j = 1; j <= 4; j++)
			{
			  face.PNum(j) = PNum(hexfaces[i - 1][j]);
			}
			break;
		}
	  }
	}

	///
	public void Invert()
	{
	  switch (GetNP())
	  {
		case 4:
		{
			netgen.GlobalMembers.Swap(ref PNum(3), ref PNum(4));
			break;
		}
		case 5:
		{
			netgen.GlobalMembers.Swap(ref PNum(1), ref PNum(4));
			netgen.GlobalMembers.Swap(ref PNum(2), ref PNum(3));
			break;
		}
		case 6:
		{
			netgen.GlobalMembers.Swap(ref PNum(1), ref PNum(4));
			netgen.GlobalMembers.Swap(ref PNum(2), ref PNum(5));
			netgen.GlobalMembers.Swap(ref PNum(3), ref PNum(6));
			break;
		}
	  }
	}

	/// split into 4 node tets
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetTets(Array<Element> & locels) const
	public void GetTets(Array<Element> locels)
	{
	  GetTetsLocal(locels);
	  int i;
	  int j;
	  for (i = 1; i <= locels.Size(); i++)
	  {
		for (j = 1; j <= 4; j++)
		{
		  locels.Elem(i).PNum(j) = PNum(locels.Elem(i).PNum(j));
		}
	  }
	}

	/// split into 4 node tets, local point nrs
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetTetsLocal(Array<Element> & locels) const
	public void GetTetsLocal(Array<Element> locels)
	{
	  int i;
	  int j;
	  locels.SetSize(0);
	  switch (GetType())
	  {
		case ELEMENT_TYPE.TET:
		{
			int[][] linels =
			{
				new int[] {1, 2, 3, 4}
			};
			for (i = 0; i < 1; i++)
			{
				Element tet = new Element(4);
				for (j = 1; j <= 4; j++)
				{
				  tet.PNum(j) = linels[i][j - 1];
				}
				locels.Append(tet);
			}
			break;
		}
		case ELEMENT_TYPE.TET10:
		{
			int[][] linels =
			{
				new int[] {1, 5, 6, 7},
				new int[] {5, 2, 8, 9},
				new int[] {6, 8, 3, 10},
				new int[] {7, 9, 10, 4},
				new int[] {5, 6, 7, 9},
				new int[] {5, 6, 9, 8},
				new int[] {6, 7, 9, 10},
				new int[] {6, 8, 10, 9}
			};
			for (i = 0; i < 8; i++)
			{
				Element tet = new Element(4);
				for (j = 1; j <= 4; j++)
				{
				  tet.PNum(j) = linels[i][j - 1];
				}
				locels.Append(tet);
			}
			break;
		}
		case ELEMENT_TYPE.PYRAMID:
		{
			int[][] linels =
			{
				new int[] {1, 2, 3, 5},
				new int[] {1, 3, 4, 5}
			};
			for (i = 0; i < 2; i++)
			{
				Element tet = new Element(4);
				for (j = 1; j <= 4; j++)
				{
				  tet.PNum(j) = linels[i][j - 1];
				}
				locels.Append(tet);
			}
			break;
		}
		case ELEMENT_TYPE.PRISM:
		case ELEMENT_TYPE.PRISM12:
		{
			int[][] linels =
			{
				new int[] {1, 2, 3, 4},
				new int[] {4, 2, 3, 5},
				new int[] {6, 5, 4, 3}
			};
			for (i = 0; i < 3; i++)
			{
				Element tet = new Element(4);
				for (j = 0; j < 4; j++)
				{
				  tet[j] = linels[i][j];
				}
				locels.Append(tet);
			}
			break;
		}
		case ELEMENT_TYPE.HEX:
		{
			int[][] linels =
			{
				new int[] {1, 7, 2, 3},
				new int[] {1, 7, 3, 4},
				new int[] {1, 7, 4, 8},
				new int[] {1, 7, 8, 5},
				new int[] {1, 7, 5, 6},
				new int[] {1, 7, 6, 2}
			};
			for (i = 0; i < 6; i++)
			{
				Element tet = new Element(4);
				for (j = 0; j < 4; j++)
				{
				  tet[j] = linels[i][j];
				}
				locels.Append(tet);
			}
			break;
		}
		default:
		{
			cerr << "GetTetsLocal not implemented for el with " << GetNP() << " nodes" << "\n";
		}
		  break;
	  }
	}

	/// returns coordinates of nodes
	// void GetNodesLocal (Array<Point<3> > & points) const;
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetNodesLocalNew(Array<Point<3>> & points) const
	public void GetNodesLocalNew(Array<Point < 3>> points)
	{
	  double[][] tetpoints =
	  {
		  new double[] {1, 0, 0},
		  new double[] {0, 1, 0},
		  new double[] {0, 0, 1},
		  new double[] {0, 0, 0}
	  };

	  double[][] prismpoints =
	  {
		  new double[] {1, 0, 0},
		  new double[] {0, 1, 0},
		  new double[] {0, 0, 0},
		  new double[] {1, 0, 1},
		  new double[] {0, 1, 1},
		  new double[] {0, 0, 1}
	  };

	  double[][] pyramidpoints =
	  {
		  new double[] {0, 0, 0},
		  new double[] {1, 0, 0},
		  new double[] {1, 1, 0},
		  new double[] {0, 1, 0},
		  new double[] {0, 0, 1}
	  };

	  double[][] tet10points =
	  {
		  new double[] {0, 0, 0},
		  new double[] {1, 0, 0},
		  new double[] {0, 1, 0},
		  new double[] {0, 0, 1},
		  new double[] {0.5, 0, 0},
		  new double[] {0, 0.5, 0},
		  new double[] {0, 0, 0.5},
		  new double[] {0.5, 0.5, 0},
		  new double[] {0.5, 0, 0.5},
		  new double[] {0, 0.5, 0.5}
	  };

	  double[][] hexpoints =
	  {
		  new double[] {0, 0, 0},
		  new double[] {1, 0, 0},
		  new double[] {1, 1, 0},
		  new double[] {0, 1, 0},
		  new double[] {0, 0, 1},
		  new double[] {1, 0, 1},
		  new double[] {1, 1, 1},
		  new double[] {0, 1, 1}
	  };



	  int np;
	  int i;
	  double[] pp = new double[3];
	  switch (GetType())
	  {
		case ELEMENT_TYPE.TET:
		{
			np = 4;
			pp = tetpoints;
			break;
		}
		case ELEMENT_TYPE.PRISM:
		case ELEMENT_TYPE.PRISM12:
		{
			np = 6;
			pp = prismpoints;
			break;
		}
		case ELEMENT_TYPE.TET10:
		{
			np = 10;
			pp = tet10points;
			break;
		}
		case ELEMENT_TYPE.PYRAMID:
		{
			np = 5;
			pp = pyramidpoints;
			break;
		}
		case ELEMENT_TYPE.HEX:
		{
			np = 8;
			pp = hexpoints;
			break;
		}
		default:
		{
			Console.Write("GetNodesLocal not impelemented for element ");
			Console.Write(GetType());
			Console.Write("\n");
			np = 0;
		pp = null;
		}
		  break;
	  }

	  points.SetSize(0);
	  for (i = 0; i < np; i++)
	  {
		points.Append(Point < 3> (pp[i][0], pp[i][1], pp[i][2]));
	  }
	}

	/// split surface into 3 node trigs
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	private int[][] GetSurfaceTriangles_tet4trigs =
	{
		new int[] {2, 3, 4},
		new int[] {3, 1, 4},
		new int[] {1, 2, 4},
		new int[] {2, 1, 3}
	};
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	private int[][] GetSurfaceTriangles_tet10trigs =
	{
		new int[] {2, 8, 9},
		new int[] {3, 10, 8},
		new int[] {4, 9, 10},
		new int[] {9, 8, 10},
		new int[] {3, 6, 10},
		new int[] {1, 7, 6},
		new int[] {4, 10, 7},
		new int[] {6, 7, 10},
		new int[] {1, 5, 7},
		new int[] {2, 9, 5},
		new int[] {4, 7, 9},
		new int[] {5, 9, 7},
		new int[] {1, 6, 5},
		new int[] {2, 5, 8},
		new int[] {3, 8, 6},
		new int[] {5, 6, 8}
	};
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	private int[][] GetSurfaceTriangles_pyramidtrigs =
	{
		new int[] {1, 3, 2},
		new int[] {1, 4, 3},
		new int[] {1, 2, 5},
		new int[] {2, 3, 5},
		new int[] {3, 4, 5},
		new int[] {4, 1, 5}
	};
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	private int[][] GetSurfaceTriangles_prismtrigs =
	{
		new int[] {1, 3, 2},
		new int[] {4, 5, 6},
		new int[] {1, 2, 4},
		new int[] {4, 2, 5},
		new int[] {2, 3, 5},
		new int[] {5, 3, 6},
		new int[] {3, 1, 6},
		new int[] {6, 1, 4}
	};
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	private int[][] GetSurfaceTriangles_hextrigs =
	{
		new int[] {1, 3, 2},
		new int[] {1, 4, 3},
		new int[] {5, 6, 7},
		new int[] {5, 7, 8},
		new int[] {1, 2, 6},
		new int[] {1, 6, 5},
		new int[] {2, 3, 7},
		new int[] {2, 7, 6},
		new int[] {3, 4, 8},
		new int[] {3, 8, 7},
		new int[] {4, 1, 8},
		new int[] {1, 5, 8}
	};

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetSurfaceTriangles(Array<Element2d> & surftrigs) const
	public void GetSurfaceTriangles(Array<Element2d> surftrigs)
	{
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	  static int tet4trigs[][3] = { { 2, 3, 4 }, { 3, 1, 4 }, { 1, 2, 4 }, { 2, 1, 3 } };

//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	  static int tet10trigs[][3] = { { 2, 8, 9 }, { 3, 10, 8}, { 4, 9, 10 }, { 9, 8, 10 }, { 3, 6, 10 }, { 1, 7, 6 }, { 4, 10, 7 }, { 6, 7, 10 }, { 1, 5, 7 }, { 2, 9, 5 }, { 4, 7, 9 }, { 5, 9, 7 }, { 1, 6, 5 }, { 2, 5, 8 }, { 3, 8, 6 }, { 5, 6, 8 } };

//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	  static int pyramidtrigs[][3] = { { 1, 3, 2 }, { 1, 4, 3 }, { 1, 2, 5 }, { 2, 3, 5 }, { 3, 4, 5 }, { 4, 1, 5 } };

//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	  static int prismtrigs[][3] = { { 1, 3, 2 }, { 4, 5, 6 }, { 1, 2, 4 }, { 4, 2, 5 }, { 2, 3, 5 }, { 5, 3, 6 }, { 3, 1, 6 }, { 6, 1, 4 } };

//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	  static int hextrigs[][3] = { { 1, 3, 2 }, { 1, 4, 3 }, { 5, 6, 7 }, { 5, 7, 8 }, { 1, 2, 6 }, { 1, 6, 5 }, { 2, 3, 7 }, { 2, 7, 6 }, { 3, 4, 8 }, { 3, 8, 7 }, { 4, 1, 8 }, { 1, 5, 8 } };

	  int j;

	  int nf;
	  int[] fp = new int[3];

	  switch (GetType())
	  {
		case ELEMENT_TYPE.TET:
		{
			nf = 4;
			fp = GetSurfaceTriangles_tet4trigs;
			break;
		}
		case ELEMENT_TYPE.PYRAMID:
		{
			nf = 6;
			fp = GetSurfaceTriangles_pyramidtrigs;
			break;
		}
		case ELEMENT_TYPE.PRISM:
		case ELEMENT_TYPE.PRISM12:
		{
			nf = 8;
			fp = GetSurfaceTriangles_prismtrigs;
			break;
		}
		case ELEMENT_TYPE.TET10:
		{
			nf = 16;
			fp = GetSurfaceTriangles_tet10trigs;
			break;
		}
		case ELEMENT_TYPE.HEX:
		{
			nf = 12;
			fp = GetSurfaceTriangles_hextrigs;
			break;
		}
		default:
		{
			nf = 0;
			fp = null;
		}
		  break;
	  }


	  surftrigs.SetSize(nf);
	  for (j = 0; j < nf; j++)
	  {
		  surftrigs.Elem(j + 1) = new Element2d(ELEMENT_TYPE.TRIG);
		  surftrigs.Elem(j + 1).PNum(1) = fp[j][0];
		  surftrigs.Elem(j + 1).PNum(2) = fp[j][1];
		  surftrigs.Elem(j + 1).PNum(3) = fp[j][2];
	  }
	}


	/// get number of 'integration points'
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetNIP() const
	public int GetNIP()
	{
	  int nip;
	  switch (typ)
	  {
		case ELEMENT_TYPE.TET:
			nip = 1;
			break;
		case ELEMENT_TYPE.TET10:
			nip = 8;
			break;
		default:
			nip = 0;
			break;
	  }
	  return nip;
	}

//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	private double[][] GetIntegrationPoint_eltetqp =
	{
		new double[] {0.25, 0.25, 0.25, 1.0 / 6.0}
	};
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	private double[][] GetIntegrationPoint_eltet10qp =
	{
		new double[] {0.585410196624969, 0.138196601125011, 0.138196601125011, 1.0 / 24.0},
		new double[] {0.138196601125011, 0.585410196624969, 0.138196601125011, 1.0 / 24.0},
		new double[] {0.138196601125011, 0.138196601125011, 0.585410196624969, 1.0 / 24.0},
		new double[] {0.138196601125011, 0.138196601125011, 0.138196601125011, 1.0 / 24.0},
		new double[] {1, 0, 0, 1},
		new double[] {0, 1, 0, 1},
		new double[] {0, 0, 1, 1},
		new double[] {0, 0, 0, 1}
	};

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetIntegrationPoint(int ip, Point<3> & p, double & weight) const
	public void GetIntegrationPoint(int ip, Point < 3> p, ref double weight)
	{
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	  static double eltetqp[1][4] = { { 0.25, 0.25, 0.25, 1.0/6.0 } };

//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	  static double eltet10qp[8][4] = { { 0.585410196624969, 0.138196601125011, 0.138196601125011, 1.0/24.0 }, { 0.138196601125011, 0.585410196624969, 0.138196601125011, 1.0/24.0 }, { 0.138196601125011, 0.138196601125011, 0.585410196624969, 1.0/24.0 }, { 0.138196601125011, 0.138196601125011, 0.138196601125011, 1.0/24.0 }, { 1, 0, 0, 1 }, { 0, 1, 0, 1 }, { 0, 0, 1, 1 }, { 0, 0, 0, 1 }};

//C++ TO C# CONVERTER TODO TASK: C# does not have an equivalent to pointers to value types:
//ORIGINAL LINE: double * pp = null;
	  double pp = null;
	  switch (typ)
	  {
		case ELEMENT_TYPE.TET:
			pp = GetIntegrationPoint_eltetqp[0][0];
			break;
		case ELEMENT_TYPE.TET10:
			pp = GetIntegrationPoint_eltet10qp[ip - 1][0];
			break;
		default:
		  throw new Exception("illegal element shape in GetIntegrationPoint");
	  }

	  p.functorMethod(0) = pp[0];
	  p.functorMethod(1) = pp[1];
	  p.functorMethod(2) = pp[2];
	  weight = pp[3];
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetTransformation(int ip, const T_POINTS & points, class DenseMatrix & trans) const;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//	void GetTransformation(int ip, T_POINTS points, DenseMatrix trans);
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetTransformation(int ip, class DenseMatrix & pmat, class DenseMatrix & trans) const
	public void GetTransformation(int ip, DenseMatrix pmat, DenseMatrix trans)
	{
	  int np = GetNP();

	  if (pmat.Width() != np || pmat.Height() != 3)
	  {
		  (*testout) << "GetTransofrmation: pmat doesn't fit" << "\n";
		  return;
	  }

	  ComputeIntegrationPointData();
	  DenseMatrix dshapep = null;
	  switch (GetType())
	  {
		case ELEMENT_TYPE.TET:
			dshapep = GlobalMembers.ipdtet.Get(ip).dshape;
			break;
		case ELEMENT_TYPE.TET10:
			dshapep = GlobalMembers.ipdtet10.Get(ip).dshape;
			break;
		default:
		  PrintSysError("Element::GetTransformation, illegal type ", (int)typ);
		  break;
	  }

	  CalcABt(pmat.functorMethod, dshapep.functorMethod, trans.functorMethod);
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetShape(const Point<3> & p, class Vector & shape) const;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//	void GetShape(Point<3> p, Vector shape);
	// void GetShapeNew (const Point<3> & p, class FlatVector & shape) const;
//C++ TO C# CONVERTER TODO TASK: The original C++ template specifier was replaced with a C# generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <typename T>
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetShapeNew(const Point<3,T> & p, TFlatVector<T> shape) const
	public void GetShapeNew<T>(Point<3,T> p, TFlatVector<T> shape)
	{
	  /*
	    if (shape.Size() < GetNP())
	    {
	    cerr << "Element::GetShape: Length not fitting" << endl;
	    return;
	    }
	  */

	  switch (typ)
	  {
		case ELEMENT_TYPE.TET:
		{
			shape.functorMethod(0) = p(0);
			shape.functorMethod(1) = p(1);
			shape.functorMethod(2) = p(2);
			shape.functorMethod(3) = 1 - p(0) - p(1) - p(2);
			break;
		}

		case ELEMENT_TYPE.TET10:
		{
			T lam1 = p(0);
			T lam2 = p(1);
			T lam3 = p(2);
			T lam4 = 1 - p(0) - p(1) - p(2);

			shape.functorMethod(0) = 2 * lam1 * (lam1 - 0.5);
			shape.functorMethod(1) = 2 * lam2 * (lam2 - 0.5);
			shape.functorMethod(2) = 2 * lam3 * (lam3 - 0.5);
			shape.functorMethod(3) = 2 * lam4 * (lam4 - 0.5);

			shape.functorMethod(4) = 4 * lam1 * lam2;
			shape.functorMethod(5) = 4 * lam1 * lam3;
			shape.functorMethod(6) = 4 * lam1 * lam4;
			shape.functorMethod(7) = 4 * lam2 * lam3;
			shape.functorMethod(8) = 4 * lam2 * lam4;
			shape.functorMethod(9) = 4 * lam3 * lam4;

			break;
		}


		case ELEMENT_TYPE.PYRAMID:
		{
			T noz = 1 - p(2);
			// if (noz == 0.0) noz = 1e-10;
			noz += T(1e-12);

			T xi = p(0) / noz;
			T eta = p(1) / noz;
			shape.functorMethod(0) = (1 - xi) * (1 - eta) * (noz);
			shape.functorMethod(1) = (xi) * (1 - eta) * (noz);
			shape.functorMethod(2) = (xi) * (eta) * (noz);
			shape.functorMethod(3) = (1 - xi) * (eta) * (noz);
			shape.functorMethod(4) = p(2);
			break;
		}
		case ELEMENT_TYPE.PYRAMID13:
		{
		T x = p(0);
		T y = p(1);
		T z = p(2);
			z *= 1 - 1e-12;
			shape.functorMethod[0] = (-z + z * (2 * x + z - 1) * (2 * y + z - 1) / (-z + 1) + (-2 * x - z + 2) * (-2 * y - z + 2)) * (-0.5 * x - 0.5 * y - 0.5 * z + 0.25);
			shape.functorMethod[1] = (0.5 * x - 0.5 * y - 0.25) * (-z - z * (2 * x + z - 1) * (2 * y + z - 1) / (-z + 1) + (2 * x + z) * (-2 * y - z + 2));
			shape.functorMethod[2] = (-z + z * (2 * x + z - 1) * (2 * y + z - 1) / (-z + 1) + (2 * x + z) * (2 * y + z)) * (0.5 * x + 0.5 * y + 0.5 * z - 0.75);
			shape.functorMethod[3] = (-0.5 * x + 0.5 * y - 0.25) * (-z - z * (2 * x + z - 1) * (2 * y + z - 1) / (-z + 1) + (2 * y + z) * (-2 * x - z + 2));
			shape.functorMethod[4] = z * (2 * z - 1);
			shape.functorMethod[5] = 2 * x * (-2 * x - 2 * z + 2) * (-2 * y - 2 * z + 2) / (-2 * z + 2);
			shape.functorMethod[6] = 4 * x * y * (-2 * x - 2 * z + 2) / (-2 * z + 2);
			shape.functorMethod[7] = 2 * y * (-2 * x - 2 * z + 2) * (-2 * y - 2 * z + 2) / (-2 * z + 2);
			shape.functorMethod[8] = 4 * x * y * (-2 * y - 2 * z + 2) / (-2 * z + 2);
			shape.functorMethod[9] = z * (-2 * x - 2 * z + 2) * (-2 * y - 2 * z + 2) / (-z + 1);
			shape.functorMethod[10] = 2 * x * z * (-2 * y - 2 * z + 2) / (-z + 1);
			shape.functorMethod[11] = 4 * x * y * z / (-z + 1);
			shape.functorMethod[12] = 2 * y * z * (-2 * x - 2 * z + 2) / (-z + 1);
			break;
		}
		case ELEMENT_TYPE.PRISM:
		{
			shape.functorMethod(0) = p(0) * (1 - p(2));
			shape.functorMethod(1) = p(1) * (1 - p(2));
			shape.functorMethod(2) = (1 - p(0) - p(1)) * (1 - p(2));
			shape.functorMethod(3) = p(0) * p(2);
			shape.functorMethod(4) = p(1) * p(2);
			shape.functorMethod(5) = (1 - p(0) - p(1)) * p(2);
			break;
		}
		case ELEMENT_TYPE.PRISM15:
		{
		T x = p(0);
		T y = p(1);
		T z = p(2);
			T lam = 1 - x - y;
			T lamz = 1 - z;
			shape.functorMethod[0] = (2 * x * x - x) * (2 * lamz * lamz - lamz);
			shape.functorMethod[1] = (2 * y * y - y) * (2 * lamz * lamz - lamz);
			shape.functorMethod[2] = (2 * lam * lam - lam) * (2 * lamz * lamz - lamz);
			shape.functorMethod[3] = (2 * x * x - x) * (2 * z * z - z);
			shape.functorMethod[4] = (2 * y * y - y) * (2 * z * z - z);
			shape.functorMethod[5] = (2 * lam * lam - lam) * (2 * z * z - z);
			shape.functorMethod[6] = 4 * x * y * (2 * lamz * lamz - lamz);
			shape.functorMethod[7] = 4 * x * lam * (2 * lamz * lamz - lamz);
			shape.functorMethod[8] = 4 * y * lam * (2 * lamz * lamz - lamz);
			shape.functorMethod[9] = x * 4 * z * (1 - z);
			shape.functorMethod[10] = y * 4 * z * (1 - z);
			shape.functorMethod[11] = lam * 4 * z * (1 - z);
			shape.functorMethod[12] = 4 * x * y * (2 * z * z - z);
			shape.functorMethod[13] = 4 * x * lam * (2 * z * z - z);
			shape.functorMethod[14] = 4 * y * lam * (2 * z * z - z);
			break;
		}
		case ELEMENT_TYPE.HEX:
		{
			shape.functorMethod(0) = (1 - p(0)) * (1 - p(1)) * (1 - p(2));
			shape.functorMethod(1) = (p(0)) * (1 - p(1)) * (1 - p(2));
			shape.functorMethod(2) = (p(0)) * (p(1)) * (1 - p(2));
			shape.functorMethod(3) = (1 - p(0)) * (p(1)) * (1 - p(2));
			shape.functorMethod(4) = (1 - p(0)) * (1 - p(1)) * (p(2));
			shape.functorMethod(5) = (p(0)) * (1 - p(1)) * (p(2));
			shape.functorMethod(6) = (p(0)) * (p(1)) * (p(2));
			shape.functorMethod(7) = (1 - p(0)) * (p(1)) * (p(2));
			break;
		}
		case ELEMENT_TYPE.HEX20:
		{
		T x = p(0);
		T y = p(1);
		T z = p(2);
		shape.functorMethod[0] = (1 - x) * (1 - y) * (1 - z);
		shape.functorMethod[1] = x * (1 - y) * (1 - z);
		shape.functorMethod[2] = x * y * (1 - z);
		shape.functorMethod[3] = (1 - x) * y * (1 - z);
		shape.functorMethod[4] = (1 - x) * (1 - y) * (z);
		shape.functorMethod[5] = x * (1 - y) * (z);
		shape.functorMethod[6] = x * y * (z);
		shape.functorMethod[7] = (1 - x) * y * (z);

			T[] sigma = {(1 - x) + (1 - y) + (1 - z), x + (1 - y) + (1 - z), x + y + (1 - z), (1 - x) + y + (1 - z), (1 - x) + (1 - y) + z, x + (1 - y) + z, x + y + z, (1 - x) + y + z};

			int[][] e =
			{
				new int[] {0, 1},
				new int[] {2, 3},
				new int[] {3, 0},
				new int[] {1, 2},
				new int[] {4, 5},
				new int[] {6, 7},
				new int[] {7, 4},
				new int[] {5, 6},
				new int[] {0, 4},
				new int[] {1, 5},
				new int[] {2, 6},
				new int[] {3, 7}
			};

			for (int i = 0; i < 12; i++)
			{
				T lame = shape.functorMethod[e[i][0]] + shape.functorMethod[e[i][1]];
				T xi = sigma[e[i][1]] - sigma[e[i][0]];
				shape.functorMethod[8 + i] = (1 - xi * xi) * lame;
			}
			for (int i = 0; i < 12; i++)
			{
				shape.functorMethod[e[i][0]] -= 0.5 * shape.functorMethod[8 + i];
				shape.functorMethod[e[i][1]] -= 0.5 * shape.functorMethod[8 + i];
			}
			break;
		}
		default:
		  throw new Exception("Element :: GetNewShape not implemented for that element");
	  }
	}

	/// matrix 2 * np
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetDShape(const Point<3> & p, class DenseMatrix & dshape) const;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//	void GetDShape(Point<3> p, DenseMatrix dshape);
//C++ TO C# CONVERTER TODO TASK: The original C++ template specifier was replaced with a C# generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <typename T>
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetDShapeNew(const Point<3,T> & p, class MatrixFixWidth<3,T> & dshape) const;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//	void GetDShapeNew<T>(Point<3,T> p, MatrixFixWidth<3,T> dshape);
	/// matrix 3 * np
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetPointMatrix(const T_POINTS & points, class DenseMatrix & pmat) const;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//	void GetPointMatrix<T>(T_POINTS points, DenseMatrix pmat);

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void ComputeIntegrationPointData() const
	public void ComputeIntegrationPointData<T>()
	{
	  switch (GetType())
	  {
		case ELEMENT_TYPE.TET:
			if (GlobalMembers.ipdtet.Size())
			{
				return;
			}
			break;
		case ELEMENT_TYPE.TET10:
			if (GlobalMembers.ipdtet10.Size())
			{
				return;
			}
			break;
		default:
		  PrintSysError("Element::ComputeIntegrationPoint, illegal type ", (int)typ);
		  break;
	  }

	  switch (GetType())
	  {
		case ELEMENT_TYPE.TET:
			GlobalMembers.ipdtet.SetSize(GetNIP());
			break;
		case ELEMENT_TYPE.TET10:
			GlobalMembers.ipdtet10.SetSize(GetNIP());
			break;
		default:
		  PrintSysError("Element::ComputeIntegrationPoint, illegal type2 ", (int)typ);
		  break;
	  }


	  for (int i = 1; i <= GetNIP(); i++)
	  {
		  IntegrationPointData ipd = new IntegrationPointData();
		  GetIntegrationPoint(i, ipd.p, ref ipd.weight);
		  ipd.shape.SetSize(GetNP());
		  ipd.dshape.SetSize(3, GetNP());

		  GetShape(ipd.p, ipd.shape);
		  GetDShape(ipd.p, ipd.dshape.functorMethod);

		  switch (GetType())
		  {
			case ELEMENT_TYPE.TET:
				GlobalMembers.ipdtet.Elem(i).reset(ipd);
				break;
			case ELEMENT_TYPE.TET10:
				GlobalMembers.ipdtet10.Elem(i).reset(ipd);
				break;
			default:
			  PrintSysError("Element::ComputeIntegrationPoint(2), illegal type ", (int)typ);
			  break;
		  }
	  }
	}


//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: double CalcJacobianBadness(const T_POINTS & points) const
	public double CalcJacobianBadness(T_POINTS points)
	{
	  int nip = GetNIP();
	  DenseMatrix trans = new DenseMatrix(3, 3);
	  DenseMatrix pmat = new DenseMatrix();

	  pmat.SetSize(3, GetNP());
	  GetPointMatrix(points, pmat.functorMethod);

	  double err = 0;
	  for (int i = 1; i <= nip; i++)
	  {
		  GetTransformation(i, pmat.functorMethod, trans.functorMethod);

		  // Frobenius norm
		  double frob = 0;
		  for (int j = 1; j <= 9; j++)
		  {
			frob += netgen.GlobalMembers.sqr(trans.Get(j));
		  }
		  frob = ngsimd.GlobalMembers.sqrt(frob);
		  frob /= 3;

		  double det = -trans.Det();

		  if (det <= 0)
		  {
			err += 1e12;
		  }
		  else
		  {
			err += frob * frob * frob / det;
		  }
	  }

	  err /= nip;
	  return err;
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: double CalcJacobianBadnessDirDeriv(const T_POINTS & points, int pi, Vec<3> & dir, double & dd) const
	public double CalcJacobianBadnessDirDeriv(T_POINTS points, int pi, Vec < 3> dir, ref double dd)
	{
	  int i;
	  int j;
	  int k;
	  int nip = GetNIP();
	  DenseMatrix trans = new DenseMatrix(3, 3);
	  DenseMatrix dtrans = new DenseMatrix(3, 3);
	  DenseMatrix hmat = new DenseMatrix(3, 3);
	  DenseMatrix pmat = new DenseMatrix();
	  DenseMatrix vmat = new DenseMatrix();

	  pmat.SetSize(3, GetNP());
	  vmat.SetSize(3, GetNP());

	  GetPointMatrix(points, pmat.functorMethod);

	  for (i = 1; i <= np; i++)
	  {
		for (j = 1; j <= 3; j++)
		{
		  vmat.Elem(j, i) = 0;
		}
	  }
	  for (j = 1; j <= 3; j++)
	  {
		vmat.Elem(j, pi) = dir.functorMethod(j - 1);
	  }



	  double err = 0;
	  dd = 0;

	  for (i = 1; i <= nip; i++)
	  {
		  GetTransformation(i, pmat.functorMethod, trans.functorMethod);
		  GetTransformation(i, vmat.functorMethod, dtrans.functorMethod);


		  // Frobenius norm
		  double frob = 0;
		  for (j = 1; j <= 9; j++)
		  {
			frob += netgen.GlobalMembers.sqr(trans.Get(j));
		  }
		  frob = ngsimd.GlobalMembers.sqrt(frob);

		  double dfrob = 0;
		  for (j = 1; j <= 9; j++)
		  {
			dfrob += trans.Get(j) * dtrans.Get(j);
		  }
		  dfrob = dfrob / frob;

		  frob /= 3;
		  dfrob /= 3;


		  double det = trans.Det();
		  double ddet = 0;

		  for (j = 1; j <= 3; j++)
		  {
			  hmat = trans.functorMethod;
			  for (k = 1; k <= 3; k++)
			  {
				hmat.Elem(k, j) = dtrans.Get(k, j);
			  }
			  ddet += hmat.Det();
		  }


		  det *= -1;
		  ddet *= -1;


		  if (det <= 0)
		  {
			err += 1e12;
		  }
		  else
		  {
			  err += frob * frob * frob / det;
			  dd += (3 * frob * frob * dfrob * det - frob * frob * frob * ddet) / (det * det);
		  }
	  }

	  err /= nip;
	  dd /= nip;
	  return err;
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: double CalcJacobianBadnessGradient(const T_POINTS & points, int pi, Vec<3> & grad) const
	public double CalcJacobianBadnessGradient(T_POINTS points, int pi, ref Vec < 3> grad)
	{
	  int nip = GetNIP();
	  DenseMatrix trans = new DenseMatrix(3, 3);
	  DenseMatrix dtrans = new DenseMatrix(3, 3);
	  DenseMatrix hmat = new DenseMatrix(3, 3);
	  DenseMatrix pmat = new DenseMatrix();
	  DenseMatrix vmat = new DenseMatrix();

	  pmat.SetSize(3, GetNP());
	  vmat.SetSize(3, GetNP());

	  GetPointMatrix(points, pmat.functorMethod);

	  for (int i = 1; i <= np; i++)
	  {
		for (int j = 1; j <= 3; j++)
		{
		  vmat.Elem(j, i) = 0;
		}
	  }
	  for (int j = 1; j <= 3; j++)
	  {
		vmat.Elem(j, pi) = 1.0;
	  }


	  double err = 0;

	  double[] dfrob = new double[3];

	  grad = 0;

	  for (int i = 1; i <= nip; i++)
	  {
		  GetTransformation(i, pmat.functorMethod, trans.functorMethod);
		  GetTransformation(i, vmat.functorMethod, dtrans.functorMethod);

		  // Frobenius norm
		  double frob = 0;
		  for (int j = 1; j <= 9; j++)
		  {
			frob += netgen.GlobalMembers.sqr(trans.Get(j));
		  }
		  frob = ngsimd.GlobalMembers.sqrt(frob);

		  for (int k = 0; k < 3; k++)
		  {
			  dfrob[k] = 0;
			  for (int j = 1; j <= 3; j++)
			  {
				dfrob[k] += trans.Get(k + 1, j) * dtrans.Get(k + 1, j);
			  }
			  dfrob[k] = dfrob[k] / (3.0 * frob);
		  }

		  frob /= 3;

		  double det = trans.Det();
		  double[] ddet = new double[3]; // = 0;

		  for (int k = 1; k <= 3; k++)
		  {
			  int km1 = (k > 1) ? (k - 1) : 3;
			  int kp1 = (k < 3) ? (k + 1) : 1;
			  ddet[k - 1] = 0;
			  for (int j = 1; j <= 3; j++)
			  {
				  int jm1 = (j > 1) ? (j - 1) : 3;
				  int jp1 = (j < 3) ? (j + 1) : 1;

				  ddet[k - 1] += (-1.0) * dtrans.Get(k, j) * (trans.Get(km1, jm1) * trans.Get(kp1, jp1) - trans.Get(km1, jp1) * trans.Get(kp1, jm1));
			  }
		  }


		  det *= -1;

		  if (det <= 0)
		  {
			err += 1e12;
		  }
		  else
		  {
			  err += frob * frob * frob / det;
			  double fac = (frob * frob) / (det * det);
			  for (int j = 0; j < 3; j++)
			  {
				grad.functorMethod(j) += fac * (3 * dfrob[j] * det - frob * ddet[j]);
			  }
		  }
	  }

	  err /= nip;
	  grad.functorMethod *= 1.0 / nip;
	  return err;
	}

	///
	// friend ostream & operator<<(ostream  & s, const Element & el);

	public void SetRefinementFlag(bool rflag = true)
	{
		flags.refflag = rflag;
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int TestRefinementFlag() const
	public int TestRefinementFlag()
	{
		return flags.refflag;
	}

	public void SetStrongRefinementFlag(bool rflag = true)
	{
		flags.strongrefflag = rflag;
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int TestStrongRefinementFlag() const
	public int TestStrongRefinementFlag()
	{
		return flags.strongrefflag;
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int Illegal() const
	public int Illegal()
	{
		return flags.illegal;
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int IllegalValid() const
	public int IllegalValid()
	{
		return flags.illegal_valid;
	}
	public void SetIllegal(int aillegal)
	{
	  flags.illegal = aillegal != 0 ? 1 : 0;
	  flags.illegal_valid = 1;
	}
	public void SetLegal(int alegal)
	{
	  flags.illegal = alegal != 0 ? 0 : 1;
	  flags.illegal_valid = 1;
	}

	public void Delete()
	{
		flags.deleted = 1;
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool IsDeleted() const
	public bool IsDeleted()
	{
#if DEBUG
	  if (pnum[0] < PointIndex.BASE && !flags.deleted)
	  {
	cerr << "Volelement has illegal pnum, but not marked as deleted" << "\n";
	  }
#endif

	  return flags.deleted;
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool IsCurved() const
	public bool IsCurved()
	{
		return is_curved;
	}
	public void SetCurved(bool acurved)
	{
		is_curved = acurved;
	}

	/*
#ifdef PARALLEL
	int GetPartition () const { return partitionNumber; }
	void SetPartition (int nr) { partitionNumber = nr; }; 
#else
	int GetPartition () const { return 0; }
#endif
	*/

	public int hp_elnr;
  }






  /**
     Edge segment.
  */
  public class Segment : System.IDisposable
  {
	///
	public Segment()
	{
		this.is_curved = false;
	  pnums[0] = -1;
	  pnums[1] = -1;
	  edgenr = -1;

	  singedge_left = 0.0;
	  singedge_right = 0.0;
	  seginfo = 0;

	  si = -1;

	  domin = -1;
	  domout = -1;
	  tlosurf = -1;

	  surfnr1 = -1;
	  surfnr2 = -1;
	  pnums[2] = -1;
	  meshdocval = 0;
	  /*
	    geominfo[0].trignum=-1;
	    geominfo[1].trignum=-1;
  
	    epgeominfo[0].edgenr = 1;
	    epgeominfo[0].dist = 0;
	    epgeominfo[1].edgenr = 1;
	    epgeominfo[1].dist = 0;
	  */

	  bcname = null;
	}

	public Segment(Segment other)
	{
		this.edgenr = other.edgenr;
		this.singedge_left = other.singedge_left;
		this.singedge_right = other.singedge_right;
		this.seginfo = other.seginfo;
		this.si = other.si;
		this.domin = other.domin;
		this.domout = other.domout;
		this.tlosurf = other.tlosurf;
		this.geominfo = new netgen.PointGeomInfo();
		this.surfnr1 = other.surfnr1;
		this.surfnr2 = other.surfnr2;
		this.epgeominfo = new netgen.EdgePointGeomInfo();
		this.meshdocval = other.meshdocval;
		this.is_curved = other.is_curved;
		this.hp_elnr = other.hp_elnr;
	  for (int j = 0; j < 3; j++)
	  {
		pnums[j] = other.pnums[j];
	  }

	  geominfo[0] = other.geominfo[0];
	  geominfo[1] = other.geominfo[1];
	  epgeominfo[0] = other.epgeominfo[0];
	  epgeominfo[1] = other.epgeominfo[1];
	  bcname = other.bcname;
	}

	public void Dispose()
	{
		;
	}

	// friend ostream & operator<<(ostream  & s, const Segment & seg);

	public PointIndex[] pnums = Arrays.InitializeWithDefaultInstances<PointIndex>(3); // p1, p2, pmid

	public int edgenr;
	///
	public double singedge_left;
	public double singedge_right;

	/// 0.. not first segment of segs, 1..first of class, 2..first of class, inverse
//C++ TO C# CONVERTER TODO TASK: C# does not allow bit fields:
	public uint seginfo:2;

	/// surface decoding index
	public int si;
	/// co dim 2 deconding index
	public int cd2i;
	/// domain number inner side
	public int domin;
	/// domain number outer side
	public int domout;
	/// top-level object number of surface
	public int tlosurf;
	///
	public PointGeomInfo[] geominfo = Arrays.InitializeWithDefaultInstances<PointGeomInfo>(2);

	/// surfaces describing edge
	public int surfnr1;
	public int surfnr2;
	///
	public EdgePointGeomInfo[] epgeominfo = Arrays.InitializeWithDefaultInstances<EdgePointGeomInfo>(2);
	///
	// int pmid; // for second order
	///
	public int meshdocval;

	// #ifdef PARALLEL
	/// number of partition for parallel computation 
	// int partitionNumber;
	// #endif

	private string bcname;
	private bool is_curved;

	/*
	  PointIndex operator[] (int i) const
	  { return (i == 0) ? p1 : p2; }

	  PointIndex & operator[] (int i) 
	  { return (i == 0) ? p1 : p2; }
	*/

//C++ TO C# CONVERTER NOTE: This 'CopyFrom' method was converted from the original copy assignment operator:
//ORIGINAL LINE: Segment& operator =(const Segment & other)
	public Segment CopyFrom(Segment other)
	{
	  if (other != this)
	  {
	  pnums[0] = other[0];
	  pnums[1] = other[1];
	  edgenr = other.edgenr;
	  singedge_left = other.singedge_left;
	  singedge_right = other.singedge_right;
	  seginfo = other.seginfo;
	  si = other.si;
	  domin = other.domin;
	  domout = other.domout;
	  tlosurf = other.tlosurf;
	  geominfo[0] = other.geominfo[0];
	  geominfo[1] = other.geominfo[1];
	  surfnr1 = other.surfnr1;
	  surfnr2 = other.surfnr2;
	  epgeominfo[0] = other.epgeominfo[0];
	  epgeominfo[1] = other.epgeominfo[1];
	  pnums[2] = other.pnums[2];
	  meshdocval = other.meshdocval;
	  hp_elnr = other.hp_elnr;
	  bcname = other.bcname;
		  is_curved = other.is_curved;
	  }

	  return this;
	}


	public int hp_elnr;

	public void SetBCName(string abcname)
	{
	  bcname = abcname;
	}

	public string BCNamePtr()
	{
		return bcname;
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const string * BCNamePtr() const
	public string BCNamePtr()
	{
		return bcname;
	}

//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	private string GetBCName_defaultstring = "default";

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const string & GetBCName() const
	public string GetBCName()
	{
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	  static string defaultstring = "default";
	  if (bcname == null)
	  {
		  return GetBCName_defaultstring;
	  }
	  return bcname;
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetNP() const
	public int GetNP()
	{
	  return (pnums[2] < 0) ? 2 : 3;
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: ELEMENT_TYPE GetType() const
	public ELEMENT_TYPE GetType()
	{
	  return (pnums[2] < 0) ? ELEMENT_TYPE.SEGMENT : ELEMENT_TYPE.SEGMENT3;
	}

	public PointIndex this [int i]
	{
		get
		{
			return new netgen.PointIndex(pnums[i]);
		}
		set
		{
			pnums[i] = value;
		}
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const PointIndex & operator [] (int i) const
	public PointIndex this [int i]
	{
		get
		{
			return new netgen.PointIndex(pnums[i]);
		}
		set
		{
			pnums[i] = value;
		}
	}


//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool IsCurved() const
	public bool IsCurved()
	{
		return is_curved;
	}
	public void SetCurved(bool acurved)
	{
		is_curved = acurved;
	}

	/*
#ifdef PARALLEL
	int GetPartition () const { return partitionNumber; }
	void SetPartition (int nr) { partitionNumber = nr; }; 
#else
	int GetPartition () const { return 0; }
#endif
	*/

	public void DoArchive(Archive ar)
	{
	  ar[] pnums[0] pnums[1] pnums[2] & edgenr & singedge_left & singedge_right & si & cd2i & domin & domout & tlosurf & surfnr1 & surfnr2 & bcname & epgeominfo[0].edgenr & epgeominfo.edgenr = Arrays.InitializeWithDefaultInstances<ar>(1);
	}
  }


  public class Element0d
  {
	public PointIndex pnum = new PointIndex();
	public string name;
	public int index;
//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//	Element0d() = default;
	public Element0d(PointIndex _pnum, int _index)
	{
		this.pnum = new netgen.PointIndex(_pnum);
		this.index = _index;
		;
	}
  }

  // class Surface;  
  // class FaceDescriptor;

  ///
  public class FaceDescriptor
  {
	/// which surface, 0 if not available
	private int surfnr;
	/// domain nr inside
	private int domin;
	/// domain nr outside
	private int domout;
	/// top level object number of surface
	private int tlosurf;
	/// boundary condition property
	private int bcprop;
	// Philippose - 06/07/2009
	// Add capability to store surface colours along with 
	// other face data
	/// surface colour (Default: R=0.0 ; G=1.0 ; B=0.0)
	private Vec3d surfcolour = new Vec3d();

	///
	private static string default_bcname = "default";
	private string bcname = default_bcname;
	/// root of linked list 
	private SurfaceElementIndex firstelement = new SurfaceElementIndex();

	private double domin_singular;
	private double domout_singular;

	public FaceDescriptor()
	{
	  surfnr = domin = domout = bcprop = 0;
	  domin_singular = domout_singular = 0.0;
	  // Philippose - 06/07/2009
	  // Initialise surface colour
	  surfcolour = new Vec3d(0.0, 1.0, 0.0);
	  tlosurf = -1;
	  // bcname = 0;
	  firstelement = -1;
	}

	public FaceDescriptor(int surfnri, int domini, int domouti, int tlosurfi)
	{
	  surfnr = surfnri;
	  domin = domini;
	  domout = domouti;
	  // Philippose - 06/07/2009
	  // Initialise surface colour
	  surfcolour = new Vec3d(0.0, 1.0, 0.0);
	  tlosurf = tlosurfi;
	  bcprop = surfnri;
	  domin_singular = domout_singular = 0.0;
	  // bcname = 0;
	  firstelement = -1;
	}

	public FaceDescriptor(Segment seg)
	{
	  surfnr = seg.si;
	  domin = seg.domin + 1;
	  domout = seg.domout + 1;
	  // Philippose - 06/07/2009
	  // Initialise surface colour
	  surfcolour = new Vec3d(0.0, 1.0, 0.0);
	  tlosurf = seg.tlosurf + 1;
	  bcprop = 0;
	  domin_singular = domout_singular = 0.0;
	  // bcname = 0;
	  firstelement = -1;
	}

	public FaceDescriptor(FaceDescriptor other)
	{
		this.surfnr = other.surfnr;
		this.domin = other.domin;
		this.domout = other.domout;
		this.tlosurf = other.tlosurf;
		this.bcprop = other.bcprop;
		this.surfcolour = new netgen.Vec3d(other.surfcolour);
		this.bcname = other.bcname;
		this.domin_singular = other.domin_singular;
		this.domout_singular = other.domout_singular;
	  firstelement = -1;
	}

	public DLL_HEADER~FaceDescriptor()
	{
		;
	}

	public int SegmentFits(Segment seg)
	{
	  return surfnr == seg.si && domin == seg.domin + 1 && domout == seg.domout + 1 && tlosurf == seg.tlosurf + 1;
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int SurfNr() const
	public int SurfNr()
	{
		return surfnr;
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int DomainIn() const
	public int DomainIn()
	{
		return domin;
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int DomainOut() const
	public int DomainOut()
	{
		return domout;
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int TLOSurface() const
	public int TLOSurface()
	{
		return tlosurf;
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int BCProperty() const
	public int BCProperty()
	{
		return bcprop;
	}


//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: double DomainInSingular() const
	public double DomainInSingular()
	{
		return domin_singular;
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: double DomainOutSingular() const
	public double DomainOutSingular()
	{
		return domout_singular;
	}

	// Philippose - 06/07/2009
	// Get Surface colour
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: Vec3d SurfColour() const
	public Vec3d SurfColour()
	{
		return new netgen.Vec3d(surfcolour);
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: DLL_HEADER const string & GetBCName() const
	public DLL_HEADER string GetBCName()
	{
		return bcname;
	}
	// string * BCNamePtr () { return bcname; }
	// const string * BCNamePtr () const  { return bcname; }
	public void SetSurfNr(int sn)
	{
		surfnr = sn;
	}
	public void SetDomainIn(int di)
	{
		domin = di;
	}
	public void SetDomainOut(int dom)
	{
		domout = dom;
	}
	public void SetBCProperty(int bc)
	{
		bcprop = bc;
	}

	/*
	const string & FaceDescriptor :: GetBCName () const
	{
	  static string defaultstring = "default";
	  if (bcname) return *bcname;
	  return defaultstring;
	}
	*/

	public void SetBCName(string bcn)
	{
	  if (bcn != null)
	  {
		bcname = bcn;
	  }
	  else
	  {
		bcn = default_bcname;
	  }
	}

	// Philippose - 06/07/2009
	// Set the surface colour
	public void SetSurfColour(Vec3d colour)
	{
		surfcolour.CopyFrom(colour);
	}

	public void SetDomainInSingular(double v)
	{
		domin_singular = v;
	}
	public void SetDomainOutSingular(double v)
	{
		domout_singular = v;
	}

	public SurfaceElementIndex FirstElement()
	{
//C++ TO C# CONVERTER TODO TASK: The following line was determined to contain a copy constructor call - this should be verified and a copy constructor should be created:
//ORIGINAL LINE: return firstelement;
		return new netgen.SurfaceElementIndex(firstelement);
	}
	// friend ostream & operator<<(ostream  & s, const FaceDescriptor & fd);
//C++ TO C# CONVERTER TODO TASK: C# has no concept of a 'friend' class:
//	friend class Mesh;

	public void DoArchive(Archive ar)
	{
	  ar & surfnr & domin & domout & tlosurf & bcprop & surfcolour.X() & surfcolour.Y() & surfcolour.Z() & bcname & domin_singular & domout_singular;
		// don't need:  firstelement
	}
  }


  public class EdgeDescriptor
  {
	private int tlosurf;
	private int[] surfnr = new int[2];
	public EdgeDescriptor()
	{
		this.tlosurf = -1;
		surfnr[0] = surfnr[1] = -1;
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int SurfNr(int i) const
	public int SurfNr(int i)
	{
		return surfnr[i];
	}
	public void SetSurfNr(int i, int nr)
	{
		surfnr[i] = nr;
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int TLOSurface() const
	public int TLOSurface()
	{
		return tlosurf;
	}
	public void SetTLOSurface(int nr)
	{
		tlosurf = nr;
	}
  }



//C++ TO C# CONVERTER WARNING: The original type declaration contained unconverted modifiers:
//ORIGINAL LINE: class DLL_HEADER MeshingParameters
  public class MeshingParameters
  {
	/**
	   3d optimization strategy:
	   // m .. move nodes
	   // M .. move nodes, cheap functional
	   // s .. swap faces
	   // c .. combine elements
	   // d .. divide elements
	   // p .. plot, no pause
	   // P .. plot, Pause
	   // h .. Histogramm, no pause
	   // H .. Histogramm, pause
	   */
	public string optimize3d = "cmdmustm";
	/// number of 3d optimization steps
	public int optsteps3d = 3;
	/**
	   2d optimization strategy:
	   // s .. swap, opt 6 lines/node
	   // S .. swap, optimal elements
	   // m .. move nodes
	   // p .. plot, no pause
	   // P .. plot, pause
	   // c .. combine
	   **/
	public string optimize2d = "smsmsmSmSmSm";
	/// number of 2d optimization steps
	public int optsteps2d = 3;
	/// power of error (to approximate max err optimization)
	public double opterrpow = 2;
	/// do block filling ?  
	public int blockfill = 1;
	/// block filling up to distance
	public double filldist = 0.1;
	/// radius of local environment (times h)
	public double safety = 5;
	/// radius of active environment (times h)
	public double relinnersafety = 3;
	/// use local h ?
	public int uselocalh = 1;
	/// grading for local h
	public double grading = 0.3;
	/// use delaunay meshing
	public int delaunay = 1;
	/// maximal mesh size
	public double maxh = 1e10;
	/// minimal mesh size
	public double minh = 0.0;
	/// file for meshsize
	public string meshsizefilename = "";
	/// start surfacemeshing from everywhere in surface
	public int startinsurface = 0;
	/// check overlapping surfaces (debug)
	public int checkoverlap = 1;
	/// check overlapping surface mesh before volume meshing
	public int checkoverlappingboundary = 1;
	/// check chart boundary (sometimes too restrictive)
	public int checkchartboundary = 1;
	/// safety factor for curvatures (elements per radius)
	public double curvaturesafety = 2;
	/// minimal number of segments per edge
	public double segmentsperedge = 1;
	/// use parallel threads
	public int parthread = 0;
	/// weight of element size w.r.t element shape
	public double elsizeweight = 0.2;
	/// init with default values

	/// start at step
	public int perfstepsstart = 0;
	/// end at step
	public int perfstepsend = 6;


	/// from mp3:
	/// give up quality class, 2d meshing
	public int giveuptol2d = 200;
	/// give up quality class, 3d meshing
	public int giveuptol = 10;
	/// maximal outer steps
	public int maxoutersteps = 10;
	/// class starting star-shape filling
	public int starshapeclass = 5;
	/// if non-zero, baseelement must have baseelnp points
	public int baseelnp = 0;
	/// quality tolerances are handled less careful
	public int sloppy = 1;

	/// limit for max element angle (150-180)
	public double badellimit = 175;

	public bool check_impossible = false;

	public int only3D_domain_nr = 0;

	///
	public int secondorder = 0;
	/// high order element curvature
	public int elementorder = 1;
	/// quad-dominated surface meshing
	public int quad = 0;
	///
	public bool try_hexes = false;
	///
	public int inverttets = 0;
	///
	public int inverttrigs = 0;
	///
	public int autozrefine = 0;
	///
	public MeshingParameters()
	{
	  // optimize3d = "cmdmustm";
	  //optimize3d = "cmdmstm";
	  // optsteps3d = 3;
	  // optimize2d = "smsmsmSmSmSm";
	  // optsteps2d = 3;
	  // opterrpow = 2;
	  // blockfill = 1;
	  // filldist = 0.1;
	  // safety = 5;
	  // relinnersafety = 3;
	  // uselocalh = 1;
	  // grading = 0.3;
	  // delaunay = 1;
	  // maxh = 1e10;
	  // minh = 0;
	  // meshsizefilename = NULL;
	  // startinsurface = 0;
	  // checkoverlap = 1;
	  // checkoverlappingboundary = 1;
	  // checkchartboundary = 1;
	  // curvaturesafety = 2;
	  // segmentsperedge = 1;
	  // parthread = 0;

	  // elsizeweight = 0.2;
	  // giveuptol2d = 200;
	  // giveuptol = 10;
	  // maxoutersteps = 10;
	  // starshapeclass = 5;
	  // baseelnp = 0;
	  // sloppy = 1;

	  // badellimit = 175;
	  // check_impossible = 0;
	  // secondorder = 0;
	}

	///
//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//	MeshingParameters(const MeshingParameters & mp2) = default;
//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//	MeshingParameters(MeshingParameters && mp2) = default;
//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//	MeshingParameters & operator = (const MeshingParameters & mp2) = default;
//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//	MeshingParameters & operator = (MeshingParameters && mp2) = default;
	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void Print(ostream & ost) const
	public void Print(ostream ost)
	{
	  ost << "Meshing parameters: " << "\n" << "optimize3d = " << optimize3d << "\n" << "optsteps3d = " << optsteps3d << "\n" << " optimize2d = " << optimize2d << "\n" << " optsteps2d = " << optsteps2d << "\n" << " opterrpow = " << opterrpow << "\n" << " blockfill = " << blockfill << "\n" << " filldist = " << filldist << "\n" << " safety = " << safety << "\n" << " relinnersafety = " << relinnersafety << "\n" << " uselocalh = " << uselocalh << "\n" << " grading = " << grading << "\n" << " delaunay = " << delaunay << "\n" << " maxh = " << maxh << "\n" << " meshsizefilename = " << meshsizefilename << "\n" << " startinsurface = " << startinsurface << "\n" << " checkoverlap = " << checkoverlap << "\n" << " checkchartboundary = " << checkchartboundary << "\n" << " curvaturesafety = " << curvaturesafety << "\n" << " segmentsperedge = " << segmentsperedge << "\n" << " parthread = " << parthread << "\n" << " elsizeweight = " << elsizeweight << "\n" << " giveuptol2d = " << giveuptol2d << "\n" << " giveuptol = " << giveuptol << "\n" << " maxoutersteps = " << maxoutersteps << "\n" << " starshapeclass = " << starshapeclass << "\n" << " baseelnp        = " << baseelnp << "\n" << " sloppy = " << sloppy << "\n" << " badellimit = " << badellimit << "\n" << " secondorder = " << secondorder << "\n" << " elementorder = " << elementorder << "\n" << " quad = " << quad << "\n" << " inverttets = " << inverttets << "\n" << " inverttrigs = " << inverttrigs << "\n";
	}

	/// 
	// void CopyFrom(const MeshingParameters & other);

	public class MeshSizePoint
	{
	  public Point < 3> pnt;
	  public double h;
	  public MeshSizePoint(Point < 3> _pnt, double _h)
	  {
		  this.pnt = _pnt.functorMethod;
		  this.h = _h;
		  ;
	  }
//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//	  MeshSizePoint() = default;
//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//	  MeshSizePoint(const MeshSizePoint &) = default;
//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//	  MeshSizePoint(MeshSizePoint &&) = default;
//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//	  MeshSizePoint & operator = (const MeshSizePoint &) = default;
//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//	  MeshSizePoint & operator = (MeshSizePoint &&) = default;
	}
	public Array<MeshSizePoint> meshsize_points = new Array<MeshSizePoint>();

	public delegate void render_functionDelegate(bool UnnamedParameter);
	public render_functionDelegate render_function;
	public void Render(bool blocking = false)
	{
	  if (render_function != null)
	  {
		render_function(blocking);
	  }
	}
  }

  public class DebugParameters
  {
	///
	public int debugoutput;
	/// use slow checks
	public int slowchecks;
	///
	public int haltsuccess;
	///
	public int haltnosuccess;
	///
	public int haltlargequalclass;
	///
	public int haltsegment;
	///
	public int haltnode;
	///
	public int haltsegmentp1;
	///
	public int haltsegmentp2;
	///
	public int haltexistingline;
	///
	public int haltoverlap;
	///
	public int haltface;
	///
	public int haltfacenr;
	///

	/*
	void MeshingParameters :: CopyFrom(const MeshingParameters & other)
	{
	  //strcpy(optimize3d,other.optimize3d); 
	  optimize3d = other.optimize3d;
	  optsteps3d = other.optsteps3d;
	  //strcpy(optimize2d,other.optimize2d); 
	  optimize2d = other.optimize2d;
	  optsteps2d = other.optsteps2d;
	  opterrpow = other.opterrpow;
	  blockfill = other.blockfill;
	  filldist = other.filldist;
	  safety = other.safety;
	  relinnersafety = other.relinnersafety;
	  uselocalh = other.uselocalh;
	  grading = other.grading;
	  delaunay = other.delaunay;
	  maxh = other.maxh;
	  //strcpy(const_cast<char*>(meshsizefilename), other.meshsizefilename);
	  //const_cast<char*>(meshsizefilename) = other.meshsizefilename; //???
	  meshsizefilename = other.meshsizefilename;
	  startinsurface = other.startinsurface;
	  checkoverlap = other.checkoverlap;
	  checkoverlappingboundary = other.checkoverlappingboundary;
	  checkchartboundary = other.checkchartboundary;
	  curvaturesafety = other.curvaturesafety;
	  segmentsperedge = other.segmentsperedge;
	  parthread = other.parthread;
	  elsizeweight = other.elsizeweight;
	  giveuptol2d = other.giveuptol2d;
	  giveuptol = other.giveuptol;
	  maxoutersteps = other.maxoutersteps;
	  starshapeclass = other.starshapeclass;
	  baseelnp = other.baseelnp;       
	  sloppy = other.sloppy;
	  badellimit = other.badellimit;
	  secondorder = other.secondorder;
	  elementorder = other.elementorder;
	  quad = other.quad;
	  inverttets = other.inverttets;
	  inverttrigs = other.inverttrigs;
	}
	*/

	public DebugParameters()
	{
	  slowchecks = 0;
	  haltsuccess = 0;
	  haltnosuccess = 0;
	  haltlargequalclass = 0;
	  haltsegment = 0;
	  haltsegmentp1 = 0;
	  haltsegmentp2 = 0;
	}
  }







  /**
     Identification of periodic surfaces, close surfaces, etc. 
  */
  public class Identifications : System.IDisposable
  {
	public enum ID_TYPE : byte
	{
		UNDEFINED = 1,
		PERIODIC = 2,
		CLOSESURFACES = 3,
		CLOSEEDGES = 4
	}


	private Mesh mesh;

	/// identify points (thin layers, periodic b.c.)  
	private INDEX_2_HASHTABLE<int> identifiedpoints = new INDEX_2_HASHTABLE<int>();

	/// the same, with info about the id-nr
	private INDEX_3_HASHTABLE<int> identifiedpoints_nr = new INDEX_3_HASHTABLE<int>();

	/// sorted by identification nr
	private TABLE<INDEX_2> idpoints_table = new TABLE<INDEX_2>();

	private Array<ID_TYPE> type = new Array<ID_TYPE>();

	/// number of identifications (or, actually used identifications ?)
	private int maxidentnr;

	///
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//	DLL_HEADER Identifications(Mesh amesh);
	///
	public void Dispose()
	{
	  ;
	  // delete identifiedpoints;
	  // delete identifiedpoints_nr;
	}

	public void Delete()
	{
	  identifiedpoints.DeleteData();
	  identifiedpoints_nr.DeleteData();

	  /*
	  delete identifiedpoints;
	  identifiedpoints = new INDEX_2_HASHTABLE<int>(100);
	  delete identifiedpoints_nr;
	  identifiedpoints_nr = new INDEX_3_HASHTABLE<int>(100);
	  */
	  maxidentnr = 0;
	}

	/*
	  Identify points pi1 and pi2, due to
	  identification nr identnr
	*/
	public void Add(PointIndex pi1, PointIndex pi2, int identnr)
	{
	  //  (*testout) << "Identification::Add, pi1 = " << pi1 << ", pi2 = " << pi2 << ", identnr = " << identnr << endl;
	  INDEX_2 pair = new INDEX_2(pi1, pi2);
	  identifiedpoints.Set(pair, identnr);

	  INDEX_3 tripl = new INDEX_3(pi1, pi2, identnr);
	  identifiedpoints_nr.Set(tripl, 1);

	  if (identnr > maxidentnr)
	  {
		  maxidentnr = identnr;
	  }

	  if (identnr + 1 > idpoints_table.Size())
	  {
		idpoints_table.ChangeSize(identnr + 1);
	  }
	  idpoints_table.Add(identnr, pair);

	  //  timestamp = NextTimeStamp();
	}


//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int Get(PointIndex pi1, PointIndex pi2) const
	public int Get(PointIndex pi1, PointIndex pi2)
	{
	  INDEX_2 pair = new INDEX_2(pi1, pi2);
	  if (identifiedpoints.Used(pair))
	  {
		return identifiedpoints.Get(pair);
	  }
	  else
	  {
		return 0;
	  }
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetSymmetric(PointIndex pi1, PointIndex pi2) const
	public int GetSymmetric(PointIndex pi1, PointIndex pi2)
	{
	  INDEX_2 pair = new INDEX_2(pi1, pi2);
	  if (identifiedpoints.Used(pair))
	  {
		return identifiedpoints.Get(pair);
	  }

	  pair = new INDEX_2(pi2, pi1);
	  if (identifiedpoints.Used(pair))
	  {
		return identifiedpoints.Get(pair);
	  }

	  return 0;
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool Get(PointIndex pi1, PointIndex pi2, int nr) const
	public bool Get(PointIndex pi1, PointIndex pi2, int nr)
	{
	  INDEX_3 tripl = new INDEX_3(pi1, pi2, nr);
	  if (identifiedpoints_nr.Used(tripl))
	  {
		return true;
	  }
	  else
	  {
		return false;
	  }
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool GetSymmetric(PointIndex pi1, PointIndex pi2, int identnr) const;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//	bool GetSymmetric(PointIndex pi1, PointIndex pi2, int identnr);

	// bool HasIdentifiedPoints() const { return identifiedpoints != nullptr; } 
	///
	public INDEX_2_HASHTABLE<int> GetIdentifiedPoints()
	{
//C++ TO C# CONVERTER TODO TASK: The following line was determined to contain a copy constructor call - this should be verified and a copy constructor should be created:
//ORIGINAL LINE: return identifiedpoints;
	  return new netgen.INDEX_2_HASHTABLE<int>(identifiedpoints);
	}

	public bool Used(PointIndex pi1, PointIndex pi2)
	{
	  return identifiedpoints.Used(new INDEX_2(pi1, pi2));
	}

	public bool UsedSymmetric(PointIndex pi1, PointIndex pi2)
	{
	  return identifiedpoints.Used(new INDEX_2(pi1, pi2)) || identifiedpoints.Used(new INDEX_2(pi2, pi1));
	}

	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetMap(int identnr, Array<int,PointIndex::BASE> & identmap, bool symmetric = false) const
	public void GetMap(int identnr, ref Array<int,PointIndex.BASE> identmap, bool symmetric = false)
	{
	  identmap.SetSize(mesh.GetNP());
	  identmap = 0;

	  if (identnr != 0)
	  {
		for (int i = 0; i < idpoints_table[identnr].Size(); i++)
		{
			INDEX_2 pair = idpoints_table[identnr][i];
			identmap[pair.I1()] = pair.I2();
			if (symmetric)
			{
			  identmap[pair.I2()] = pair.I1();
			}
		}
	  }

	  else
	  {
		  Console.Write("getmap, identnr = ");
		  Console.Write(identnr);
		  Console.Write("\n");

		  for (int i = 1; i <= identifiedpoints_nr.GetNBags(); i++)
		  {
			for (int j = 1; j <= identifiedpoints_nr.GetBagSize(i); j++)
			{
				INDEX_3 i3 = new INDEX_3();
				int dummy;
				identifiedpoints_nr.GetData(i, j, i3, dummy);

				if (i3.I3() == identnr || identnr == 0)
				{
					identmap.Elem(i3.I1()) = i3.I2();
					if (symmetric)
					{
					  identmap.Elem(i3.I2()) = i3.I1();
					}
				}
			}
		  }
	  }

	}

	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: ID_TYPE GetType(int identnr) const
	public ID_TYPE GetType(int identnr)
	{
	  if (identnr <= type.Size())
	  {
	return type[identnr - 1];
	  }
	  else
	  {
	return ID_TYPE.UNDEFINED;
	  }
	}
	public void SetType(int identnr, ID_TYPE t)
	{
	  while (type.Size() < identnr)
	  {
	type.Append(ID_TYPE.UNDEFINED);
	  }
	  type[identnr - 1] = t;
	}

	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetPairs(int identnr, Array<INDEX_2> & identpairs) const
	public void GetPairs(int identnr, Array<INDEX_2> identpairs)
	{
	  identpairs.SetSize(0);

	  if (identnr == 0)
	  {
		for (int i = 1; i <= identifiedpoints.GetNBags(); i++)
		{
		  for (int j = 1; j <= identifiedpoints.GetBagSize(i); j++)
		  {
			  INDEX_2 i2 = new INDEX_2();
			  int nr;
			  identifiedpoints.GetData(i, j, ref i2, ref nr);
			  identpairs.Append(i2);
		  }
		}
	  }
	  else
	  {
		for (int i = 1; i <= identifiedpoints_nr.GetNBags(); i++)
		{
		  for (int j = 1; j <= identifiedpoints_nr.GetBagSize(i); j++)
		  {
			  INDEX_3 i3 = new INDEX_3();
			  int dummy;
			  identifiedpoints_nr.GetData(i, j, i3, dummy);

			  if (i3.I3() == identnr)
			  {
				identpairs.Append(new INDEX_2(i3.I1(), i3.I2()));
			  }
		  }
		}
	  }
	}

	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetMaxNr() const
	public int GetMaxNr()
	{
		return maxidentnr;
	}

	/// remove secondorder
	public void SetMaxPointNr(int maxpnum)
	{
	  for (int i = 1; i <= identifiedpoints.GetNBags(); i++)
	  {
		for (int j = 1; j <= identifiedpoints.GetBagSize(i); j++)
		{
			INDEX_2 i2 = new INDEX_2();
			int nr;
			identifiedpoints.GetData(i, j, ref i2, ref nr);

			if (i2.I1() > maxpnum || i2.I2() > maxpnum)
			{
				i2.I1() = i2.I2() = -1;
				identifiedpoints.SetData(i, j, i2, -1);
			}
		}
	  }
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void Print(ostream & ost) const
	public void Print(ostream ost)
	{
	  ost << "Identifications:" << "\n";
	  ost << "pairs: " << "\n" << identifiedpoints << "\n";
	  ost << "pairs and nr: " << "\n" << identifiedpoints_nr << "\n";
	  ost << "table: " << "\n" << idpoints_table << "\n";
	}

	public void DoArchive(Archive ar)
	{
	  ar maxidentnr;
	  ar & identifiedpoints & identifiedpoints_nr;

	  ar idpoints_table;
	  if (ar.Output())
	  {
		  uint s = type.Size();
		  ar s;
		  foreach (var t in type)
		  {
			ar & (byte)(t);
		  }
	  }
	  else
	  {
		  uint s;
		  ar s;
		  type.SetSize(s);
		  foreach (var t in type)
		  {
			ar & (byte)(t);
		  }
	  }
	}
  }
}







namespace netgen
{

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void Element2d::GetTransformation(int ip, const Array<Point2d> & points, DenseMatrix & trans) const


//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void Element2d::GetShape(const Point2d & p, Vector & shape) const



//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void Element2d::GetShapeNew(const Point<2> & p, FlatVector & shape) const









//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void Element2d::GetDShape(const Point2d & p, DenseMatrix & dshape) const


//C++ TO C# CONVERTER TODO TASK: The original C++ template specifier was replaced with a C# generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <typename T>
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void Element2d::GetDShapeNew(const Point<2,T> & p, MatrixFixWidth<2,T> & dshape) const





//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void Element2d::GetPointMatrix(const Array<Point2d> & points, DenseMatrix & pmat) const


#if OLD
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void Element::GetNodesLocal(Array<Point3d> & points) const
#endif

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void Element::GetTransformation(int ip, const T_POINTS & points, DenseMatrix & trans) const


//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void Element::GetShape(const Point<3> & hp, Vector & shape) const



//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void Element::GetDShape(const Point<3> & hp, DenseMatrix & dshape) const

//C++ TO C# CONVERTER TODO TASK: The original C++ template specifier was replaced with a C# generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <typename T>
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void Element::GetDShapeNew(const Point<3,T> & p, MatrixFixWidth<3,T> & dshape) const


//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void Element::GetPointMatrix(const T_POINTS & points, DenseMatrix & pmat) const





}

