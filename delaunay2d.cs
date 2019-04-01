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

// not yet working ....

namespace netgen
{



  public class DelaunayTrig
  {
	private PointIndex[] pnums = Arrays.InitializeWithDefaultInstances<PointIndex>(3);
	private Point < 3> c;
	private double r;
	private double rad2;
	public DelaunayTrig()
	{
		;
	}
	public DelaunayTrig(int p1, int p2, int p3)
	{
		pnums[0] = p1;
		pnums[1] = p2;
		pnums[2] = p3;
	}

	public PointIndex this [int j]
	{
		get
		{
			return new netgen.PointIndex(pnums[j]);
		}
		set
		{
			pnums[j] = value;
		}
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const PointIndex & operator [] (int j) const
	public PointIndex this [int j]
	{
		get
		{
			return new netgen.PointIndex(pnums[j]);
		}
		set
		{
			pnums[j] = value;
		}
	}

	public void CalcCenter(Mesh mesh)
	{
	  Point < 3> p1 = mesh[pnums[0]];
	  Point < 3> p2 = mesh[pnums[1]];
	  Point < 3> p3 = mesh[pnums[2]];
	  Vec < 3> v1 = p2 - p1;
	  Vec < 3> v2 = p3 - p1;
	  Mat < 2,2> mat, inv;
	  mat(0,0) = v1 * v1;
	  mat(0,1) = v1 * v2;
	  mat(1,0) = v2 * v1;
	  mat(1,1) = v2 * v2;
	  Vec < 2> rhs, sol;
	  rhs(0) = 0.5 * v1 * v1;
	  rhs(1) = 0.5 * v2 * v2;
	  netgen.GlobalMembers.CalcInverse(mat, inv);
	  sol = inv * rhs;

	  c = p1 + sol(0) * v1 + sol(1) * v2;
	  rad2 = netgen.GlobalMembers.Dist2(c, p1);
	  r = ngsimd.GlobalMembers.sqrt(rad2);
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: Point<3> Center() const
	public Point < 3> Center()
	{
		return c;
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: double Radius2() const
	public double Radius2()
	{
		return rad2;
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: Box<3> BoundingBox() const
	public Box < 3> BoundingBox()
	{
		return Box < 3> (c - Vec < 3>(r,r,0.1), c + Vec < 3>(r,r,0.1));
	}
  }



}
