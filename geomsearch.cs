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
/* File:   geomsearch.hh                                                  */
/* Author: Johannes Gerstmayr                                             */
/* Date:   19. Nov. 97                                                    */
/**************************************************************************/

//C++ TO C# CONVERTER NOTE: C# has no need of forward class declarations:
//class FrontPoint3;
//C++ TO C# CONVERTER NOTE: C# has no need of forward class declarations:
//class FrontFace;
//C++ TO C# CONVERTER NOTE: C# has no need of forward class declarations:
//class MiniElement2d;

  /// class for quick access of 3D-elements; class cannot delete elements, but only append
public class GeomSearch3d
{

  ///
  public GeomSearch3d()
  {
	size.i1 = 0;
	size.i2 = 0;
	size.i3 = 0;
  }

  ///
  public virtual void Dispose()
  {
	//delete old Hashtable:
	if (size.i1 != 0)
	{
	for (int i = 0; i < size.i1 * size.i2 * size.i3; i++)
	{
	  if (hashtable[i] != null)
	  {
		  hashtable[i].Dispose();
	  }
	}
	}
  }

  ///
  public void Init(Array <FrontPoint3,PointIndex.BASE, PointIndex> pointsi, Array <FrontFace> facesi)
  {
	points = pointsi;
	faces = facesi;
	size.i1 = 0;
	size.i2 = 0;
	size.i3 = 0;
	reset = 1;
	hashcount = 1;
  }

  ///get elements max extension
  public void ElemMaxExt(Point3d minp, Point3d maxp, MiniElement2d elem)
  {
	maxp.X() = (*points)[elem.PNum(1)].P()(0);
	maxp.Y() = (*points)[elem.PNum(1)].P()(1);
	maxp.Z() = (*points)[elem.PNum(1)].P()(2);
	minp.X() = (*points)[elem.PNum(1)].P()(0);
	minp.Y() = (*points)[elem.PNum(1)].P()(1);
	minp.Z() = (*points)[elem.PNum(1)].P()(2);

	for (int i = 2; i <= 3; i++)
	{
	maxp.X() = netgen.GlobalMembers.max2((*points)[elem.PNum(i)].P()(0), maxp.X());
	maxp.Y() = netgen.GlobalMembers.max2((*points)[elem.PNum(i)].P()(1), maxp.Y());
	maxp.Z() = netgen.GlobalMembers.max2((*points)[elem.PNum(i)].P()(2), maxp.Z());
	minp.X() = netgen.GlobalMembers.min2((*points)[elem.PNum(i)].P()(0), minp.X());
	minp.Y() = netgen.GlobalMembers.min2((*points)[elem.PNum(i)].P()(1), minp.Y());
	minp.Z() = netgen.GlobalMembers.min2((*points)[elem.PNum(i)].P()(2), minp.Z());
	}
  }

  ///get minimum coordinates of two points ->p2
  public void MinCoords(Point3d p1, Point3d p2)
  {
	p2.X() = netgen.GlobalMembers.min2(p1.X(), p2.X());
	p2.Y() = netgen.GlobalMembers.min2(p1.Y(), p2.Y());
	p2.Z() = netgen.GlobalMembers.min2(p1.Z(), p2.Z());
  }

  ///get minimum coordinates of two points ->p2
  public void MaxCoords(Point3d p1, Point3d p2)
  {
	p2.X() = netgen.GlobalMembers.max2(p1.X(), p2.X());
	p2.Y() = netgen.GlobalMembers.max2(p1.Y(), p2.Y());
	p2.Z() = netgen.GlobalMembers.max2(p1.Z(), p2.Z());
  }

  ///create a hashtable from an existing array of triangles
  ///sizei = number of pieces in one direction
  public void Create()
  {
	int i;
	int j;
	int k;
	if (reset)
	{
	const double hashelemsizefactor = 4;
	reset = 0;
	/*
	  minext=Point3d(MAXDOUBLE, MAXDOUBLE, MAXDOUBLE);
	  maxext=Point3d(MINDOUBLE, MINDOUBLE, MINDOUBLE);
	*/
	ElemMaxExt(minext, maxext, faces.Get(1).Face());
	Point3d maxp = new Point3d();
	Point3d minp = new Point3d();
	Vec3d midext = new Vec3d(0, 0, 0);

	//get max Extension of Frontfaces
	for (i = 1; i <= faces.Size(); i++)
	{
		ElemMaxExt(minp, maxp, faces.Get(i).Face());
		MinCoords(minp, minext);
		MaxCoords(maxp, maxext);
		midext += maxp - minp;
	}


	maxextreal = maxext;
	maxext = maxext + 1e-4 * (maxext - minext);

	midext *= 1.0 / faces.Size();
	Vec3d boxext = maxext - minext;

	//delete old Hashtable:
	if (size.i1 != 0)
	{
		for (i = 1; i <= size.i1 * size.i2 * size.i3; i++)
		{
		hashtable.Get(i) = null;
		}
	}

	size.i1 = (int)(boxext.X() / midext.X() / hashelemsizefactor + 1);
	size.i2 = (int)(boxext.Y() / midext.Y() / hashelemsizefactor + 1);
	size.i3 = (int)(boxext.Z() / midext.Z() / hashelemsizefactor + 1);
	// PrintMessage (5, "hashsizes = ", size.i1, ", ", size.i2, ", ", size.i3);

	elemsize.X() = boxext.X() / size.i1;
	elemsize.Y() = boxext.Y() / size.i2;
	elemsize.Z() = boxext.Z() / size.i3;

	//create Hasharrays:
	hashtable.SetSize(size.i1 * size.i2 * size.i3);
	for (i = 1; i <= size.i1; i++)
	{
		for (j = 1; j <= size.i2; j++)
		{
		for (k = 1; k <= size.i3; k++)
		{
			int ind = i + (j - 1) * size.i1 + (k - 1) * size.i2 * size.i1;
			hashtable.Elem(ind) = new Array <int> ();
		}
		}
	}
	}
	else
	{
	//Clear all Hash-Arrays
	for (i = 1; i <= size.i1; i++)
	{
		for (j = 1; j <= size.i2; j++)
		{
		for (k = 1; k <= size.i3; k++)
		{
			int ind = i + (j - 1) * size.i1 + (k - 1) * size.i2 * size.i1;
			hashtable.Elem(ind).SetSize(0);
		}
		}
	}
	}

	//Faces in Hashtable einfuegen:
	for (i = 1; i <= faces.Size(); i++)
	{
	AddElem(faces.Get(i).Face(),i);
	}

  }

  ///add new element to Hashtable
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  void AddElem(MiniElement2d elem, INDEX elemnum);

  ///GetLocal faces in sphere with radius xh and middlepoint p
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  void GetLocals(Array<MiniElement2d> locfaces, Array<INDEX> findex, INDEX fstind, Point3d p0, double xh);


  private Array<FrontFace>[] faces; // Pointers to Arrays in Adfront
  private Array<FrontPoint3,PointIndex.BASE, PointIndex>[] points;

  private Array<Array <int>> hashtable = new Array<Array <int>>();

  private Point3d minext = new Point3d(); //extension of Hashdomain
  private Point3d maxext = new Point3d();
  private Point3d maxextreal = new Point3d();
  private Vec3d elemsize = new Vec3d(); //size of one Hash-Element

  private threeint size = new threeint(); // size of Hashtable in each direction
  private int reset;
  private int hashcount;
}

























































