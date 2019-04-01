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

namespace netgen
{


  public class Opti2dLocalData
  {
	public readonly MeshOptimize2d meshthis;
	public MeshPoint sp1 = new MeshPoint();
	public PointGeomInfo gi1 = new PointGeomInfo();
	public Vec < 3> normal;
	public Vec t1 = new Vec();
	public Vec t2 = new Vec();
	public Array<SurfaceElementIndex> locelements = new Array<SurfaceElementIndex>();
	public Array<int> locrots = new Array<int>();
	public Array<double> lochs = new Array<double>();
	public Array<Point < 3>> loc_pnts2 = new Array<Point < 3>>();
	public Array<Point < 3>> loc_pnts3 = new Array<Point < 3>>();
  // static int locerr2;
	public double locmetricweight;
	public double loch;
	public int surfi;
	public int surfi2;
	public int uselocalh;
	public Opti2dLocalData()
	{
	  locmetricweight = 0;
	}
  }


  public class Opti2SurfaceMinFunction : MinFunction
  {
	private readonly Mesh mesh;
	private Opti2dLocalData ld;
	public Opti2SurfaceMinFunction(Mesh amesh, Opti2dLocalData ald)
	{
		this.mesh = new netgen.Mesh(amesh);
		this.ld = new netgen.Opti2dLocalData(ald);
	}


//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual double Func(const Vector & x) const
	public override double Func(Vector x)
	{
	  Vec < 3> n;

	  double badness = 0;

	  ld.meshthis.GetNormalVector(ld.surfi, ld.sp1, ld.gi1, n);
	  Point < 3> pp1 = ld.sp1 + x(0) * ld.t1.functorMethod + x(1) * ld.t2.functorMethod;

	  for (int j = 0; j < ld.locelements.Size(); j++)
	  {
		  Vec < 3> e1 = ld.loc_pnts2[j] - pp1;
		  Vec < 3> e2 = ld.loc_pnts3[j] - pp1;

		  if (ld.uselocalh != 0)
		  {
			  ld.loch = ld.lochs[j];
		  }

		  if (netgen.GlobalMembers.Determinant(e1, e2, n) > 1e-8 * ld.loch * ld.loch)
		  {
			  badness += netgen.GlobalMembers.CalcTriangleBadness(pp1, ld.loc_pnts2[j], ld.loc_pnts3[j], ld.locmetricweight, ld.loch);
		  }
		  else
		  {
			  badness += 1e8;
		  }
	  }

	  return badness;
	}


//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual double FuncGrad(const Vector & x, Vector & g) const
	public override double FuncGrad(Vector x, Vector g)
	{
	  Vec < 3> vgrad;
	  Point < 3> pp1;

	  vgrad = 0;
	  double badness = 0;

	  pp1 = ld.sp1 + x(0) * ld.t1.functorMethod + x(1) * ld.t2.functorMethod;

	  for (int j = 0; j < ld.locelements.Size(); j++)
	  {
		  Vec < 3> e1 = ld.loc_pnts2[j] - pp1;
		  Vec < 3> e2 = ld.loc_pnts3[j] - pp1;

		  if (ld.uselocalh != 0)
		  {
			  ld.loch = ld.lochs[j];
		  }

		  if (netgen.GlobalMembers.Determinant(e1, e2, ld.normal) > 1e-8 * ld.loch * ld.loch)
		  {
			  Vec < 3> hgrad;
			  badness += netgen.GlobalMembers.CalcTriangleBadnessGrad(pp1, ld.loc_pnts2[j], ld.loc_pnts3[j], hgrad, ld.locmetricweight, ld.loch);
			  vgrad += hgrad;
		  }
		  else
		  {
			  badness += 1e8;
		  }
	  }
	  g(0) = ld.t1.functorMethod * vgrad;
	  g(1) = ld.t2.functorMethod * vgrad;
	  return badness;
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual double FuncDeriv(const Vector & x, const Vector & dir, double & deriv) const
	public override double FuncDeriv(Vector x, Vector dir, ref double deriv)
	{
	  deriv = 0;
	  double badness = 0;

	  Point < 3> pp1 = ld.sp1 + x(0) * ld.t1.functorMethod + x(1) * ld.t2.functorMethod;
	  Vec < 3> dir3d = dir(0) * ld.t1.functorMethod + dir(1) * ld.t2.functorMethod;

	  for (int j = 0; j < ld.locelements.Size(); j++)
	  {
		  Vec < 3> e1 = ld.loc_pnts2[j] - pp1;
		  Vec < 3> e2 = ld.loc_pnts3[j] - pp1;

		  if (ld.uselocalh != 0)
		  {
			  ld.loch = ld.lochs[j];
		  }

		  if (netgen.GlobalMembers.Determinant(e1, e2, ld.normal) > 1e-8 * ld.loch * ld.loch)
		  {
			  Vec < 3> hgrad;
			  badness += netgen.GlobalMembers.CalcTriangleBadnessGrad(pp1, ld.loc_pnts2[j], ld.loc_pnts3[j], hgrad, ld.locmetricweight, ld.loch);
			  deriv += dir3d * hgrad;
		  }
		  else
		  {
			  badness += 1e8;
		  }
	  }

	  // cout << "deriv = " << deriv << " =?= ";
	  return badness;
	  /*
	  static int timer = NgProfiler::CreateTimer ("opti2surface - deriv");
	  NgProfiler::RegionTimer reg (timer);

	  double eps = 1e-6;
	  Vector xr(2), xl(2);
	  xr = x; xl = x;
	  for (int i = 0; i < 2; i++)
	    {
	      xr(i) = x(i) + eps * dir(i);
	      xl(i) = x(i) - eps * dir(i);
	    }
	  deriv = (Func (xr) - Func(xl) ) / (2*eps); 
	  cout << deriv << endl;
	  return Func(x);
	  */
	}




	/*
	double Opti2SurfaceMinFunction :: 
	Func (const Vector & x) const
	{
	  static int timer = NgProfiler::CreateTimer ("opti2surface - func");
	  NgProfiler::RegionTimer reg (timer);
  
	  Vector g(x.Size());
	  return FuncGrad (x, g);
	}
	*/

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual double XXFuncGrad(const Vector & x, Vector & grad) const
	public virtual double XXFuncGrad(Vector x, Vector grad)
	{
	  // static int timer = NgProfiler::CreateTimer ("opti2surface - funcgrad");
	  // NgProfiler::RegionTimer reg (timer);

	  Vec < 3> n, vgrad;
	  Point < 3> pp1;

	  vgrad = 0;
	  double badness = 0;

	  ld.meshthis.GetNormalVector(ld.surfi, ld.sp1, ld.gi1, n);
	  pp1 = ld.sp1 + x(0) * ld.t1.functorMethod + x(1) * ld.t2.functorMethod;

	  //  meshthis -> ProjectPoint (surfi, pp1);
	  // meshthis -> GetNormalVector (surfi, pp1, n);

	  for (int j = 0; j < ld.locelements.Size(); j++)
	  {
		  double g1x;
		  double g1y;
		  double hbadness;

		  Vec < 3> e1 = ld.loc_pnts2[j] - pp1;
		  Vec < 3> e2 = ld.loc_pnts3[j] - pp1;

	  if (ld.uselocalh != 0)
	  {
		  ld.loch = ld.lochs[j];
	  }

	  double e1 = e1.Length();
	  if (netgen.GlobalMembers.Determinant(e1, e2, n) > 1e-8 * e1 * e2.Length())
	  {
		  e1 /= e1;
		  double e1e2 = e1 * e2;
			  e2 -= e1e2 * e1;
		  double e2 = e2.Length();

			  netgen.GlobalMembers.CalcTriangleBadness(e1, e1e2, e2, ld.locmetricweight, ld.loch, hbadness, g1x, g1y);

		  badness += hbadness;
			  vgrad += g1x * e1 + (g1y / e2) * e2;
	  }
	  else
	  {
		  // (*testout) << "very very bad badness" << endl;
		  badness += 1e8;
	  }
	  }

	  // vgrad -=  (vgrad * n) * n;
	  grad(0) = vgrad * ld.t1.functorMethod;
	  grad(1) = vgrad * ld.t2.functorMethod;
	  return badness;
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual double XXFuncDeriv(const Vector & x, const Vector & dir, double & deriv) const
	public virtual double XXFuncDeriv(Vector x, Vector dir, ref double deriv)
	{
	  // static int timer = NgProfiler::CreateTimer ("opti2surface - funcderiv");
	  // NgProfiler::RegionTimer reg (timer);

	  Vec < 3> n, vgrad;
	  Point < 3> pp1;

	  vgrad = 0;
	  double badness = 0;

	  ld.meshthis.GetNormalVector(ld.surfi, ld.sp1, ld.gi1, n);
	  pp1 = ld.sp1 + x(0) * ld.t1.functorMethod + x(1) * ld.t2.functorMethod;

	  for (int j = 0; j < ld.locelements.Size(); j++)
	  {
		  double g1x;
		  double g1y;
		  double hbadness;

		  /*
		  int roti = ld.locrots[j];
		  const Element2d & bel = mesh[ld.locelements[j]];
	  Vec<3> e1 = mesh[bel.PNumMod(roti + 1)] - pp1;
	  Vec<3> e2 = mesh[bel.PNumMod(roti + 2)] - pp1;
		  */
		  Vec < 3> e1 = ld.loc_pnts2[j] - pp1;
		  Vec < 3> e2 = ld.loc_pnts3[j] - pp1;
	  if (ld.uselocalh != 0)
	  {
		  ld.loch = ld.lochs[j];
	  }

	  double e1 = e1.Length();
	  if (netgen.GlobalMembers.Determinant(e1, e2, n) > 1e-8 * e1 * e2.Length())
	  {
		  e1 /= e1;
		  double e1e2 = e1 * e2;
		  e2 -= e1e2 * e1;
		  double e2 = e2.Length();
		  netgen.GlobalMembers.CalcTriangleBadness(e1, e1e2, e2, ld.locmetricweight, ld.loch, hbadness, g1x, g1y);

		  badness += hbadness;
			  vgrad += g1x * e1 + (g1y / e2) * e2;
	  }
	  else
	  {
		  // (*testout) << "very very bad badness" << endl;
		  badness += 1e8;
	  }
	  }

	  // vgrad -= (vgrad * n) * n;
	  deriv = dir(0) * (vgrad * ld.t1.functorMethod) + dir(1) * (vgrad * ld.t2.functorMethod);

	  return badness;
	}

  }












  public class Opti2EdgeMinFunction : MinFunction
  {
	private readonly Mesh mesh;
	private Opti2dLocalData ld;

	public Opti2EdgeMinFunction(Mesh amesh, Opti2dLocalData ald)
	{
		this.mesh = new netgen.Mesh(amesh);
		this.ld = new netgen.Opti2dLocalData(ald);
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual double FuncGrad(const Vector & x, Vector & grad) const
	public override double FuncGrad(Vector x, Vector grad)
	{
	  int j;
	  int rot;
	  Vec < 3> n1, n2, v1, v2, e1, e2, vgrad;
	  Point < 3> pp1;
	  Vec < 2> g1;
	  double badness;
	  double hbadness;

	  vgrad = 0.0;
	  badness = 0;

	  pp1 = ld.sp1 + x(0) * ld.t1.functorMethod;
	  ld.meshthis.ProjectPoint2(ld.surfi, ld.surfi2, ref pp1);

	  for (j = 0; j < ld.locelements.Size(); j++)
	  {
	  rot = ld.locrots[j];
	  Element2d bel = mesh[ld.locelements[j]];

	  v1 = mesh[bel.PNumMod(rot + 1)] - pp1;
	  v2 = mesh[bel.PNumMod(rot + 2)] - pp1;

	  e1 = v1;
	  e2 = v2;
	  e1 /= e1.Length();
	  e2 -= (e1 * e2) * e1;
	  e2 /= e2.Length();

	  if (ld.uselocalh != 0)
	  {
		  ld.loch = ld.lochs[j];
	  }
	  netgen.GlobalMembers.CalcTriangleBadness((e1 * v1), (e1 * v2), (e2 * v2), ld.locmetricweight, ld.loch, hbadness, g1(0), g1(1));

	  badness += hbadness;
		  vgrad += g1(0) * e1 + g1(1) * e2;
	  }

	  ld.meshthis.GetNormalVector(ld.surfi, pp1, n1);
	  ld.meshthis.GetNormalVector(ld.surfi2, pp1, n2);

	  v1 = netgen.GlobalMembers.Cross(n1, n2);
	  v1.Normalize();

	  grad(0) = (vgrad * v1) * (ld.t1.functorMethod * v1);

	  return badness;
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual double Func(const Vector & x) const
	public override double Func(Vector x)
	{
	  Vector g = new Vector(x.Size());
	  return FuncGrad(x, g);
	}
  }




  public class Opti2SurfaceMinFunctionJacobian : MinFunction
  {
	private readonly Mesh mesh;
	private Opti2dLocalData ld;

	public Opti2SurfaceMinFunctionJacobian(Mesh amesh, Opti2dLocalData ald)
	{
		this.mesh = new netgen.Mesh(amesh);
		this.ld = new netgen.Opti2dLocalData(ald);
	}
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	private Array<Point2d> FuncGrad_pts2d = new Array<Point2d>();

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual double FuncGrad(const Vector & x, Vector & grad) const
	public override double FuncGrad(Vector x, ref Vector grad)
	{
	  // from 2d:

	  int lpi;
	  int gpi;
	  Vec < 3> n, vgrad;
	  Point < 3> pp1;
	  Vec2d g1 = new Vec2d();
	  Vec2d vdir = new Vec2d();
	  double badness;
	  double hbad;
	  double hderiv;

	  vgrad = 0;
	  badness = 0;

	  ld.meshthis.GetNormalVector(ld.surfi, ld.sp1, ld.gi1, n);

	  pp1 = ld.sp1 + x(0) * ld.t1.functorMethod + x(1) * ld.t2.functorMethod;

	  //  meshthis -> ProjectPoint (surfi, pp1);
	  //  meshthis -> GetNormalVector (surfi, pp1, n);

//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	  static Array<Point2d> pts2d;
	  FuncGrad_pts2d.SetSize(mesh.GetNP());

	  grad = 0;

	  for (int j = 1; j <= ld.locelements.Size(); j++)
	  {
	  lpi = ld.locrots.Get(j);
	  Element2d bel = mesh[ld.locelements.Get(j)];

	  gpi = bel.PNum(lpi);

	  for (int k = 1; k <= bel.GetNP(); k++)
	  {
		  PointIndex pi = bel.PNum(k);
		  FuncGrad_pts2d.Elem(pi) = new Point2d(ld.t1.functorMethod * (new mesh.Point(pi) - ld.sp1), ld.t2.functorMethod * (new mesh.Point(pi) - ld.sp1));
	  }
	  FuncGrad_pts2d.Elem(gpi) = new Point2d(x(0), x(1));


	  for (int k = 1; k <= 2; k++)
	  {
		  if (k == 1)
		  {
			vdir = new Vec2d(1, 0);
		  }
		  else
		  {
			vdir = new Vec2d(0, 1);
		  }

		  hbad = bel.CalcJacobianBadnessDirDeriv(FuncGrad_pts2d, lpi, vdir, ref hderiv);

		  grad(k - 1) += hderiv;
		  if (k == 1)
		  {
			badness += hbad;
		  }
	  }
	  }


	  /*
	    vgrad.Add (-(vgrad * n), n);
  
	    grad.Elem(1) = vgrad * t1;
	    grad.Elem(2) = vgrad * t2;
	  */
	  return badness;
	}

//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	private Array<Point2d> FuncDeriv_pts2d = new Array<Point2d>();

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual double FuncDeriv(const Vector & x, const Vector & dir, double & deriv) const
	public override double FuncDeriv(Vector x, Vector dir, ref double deriv)
	{
	  // from 2d:

	  int j;
	  int k;
	  int lpi;
	  int gpi;
	  Vec < 3> n, vgrad;
	  Point < 3> pp1;
	  Vec2d g1 = new Vec2d();
	  Vec2d vdir = new Vec2d();
	  double badness;
	  double hbad;
	  double hderiv;

	  vgrad = 0;
	  badness = 0;

	  ld.meshthis.GetNormalVector(ld.surfi, ld.sp1, ld.gi1, n);

	  // pp1 = sp1;
	  //    pp1.Add2 (x.Get(1), t1, x.Get(2), t2);
	  pp1 = ld.sp1 + x(0) * ld.t1.functorMethod + x(1) * ld.t2.functorMethod;

//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	  static Array<Point2d> pts2d;
	  FuncDeriv_pts2d.SetSize(mesh.GetNP());

	  deriv = 0;

	  for (j = 1; j <= ld.locelements.Size(); j++)
	  {
	  lpi = ld.locrots.Get(j);
	  Element2d bel = mesh[ld.locelements.Get(j)];

	  gpi = bel.PNum(lpi);

	  for (k = 1; k <= bel.GetNP(); k++)
	  {
		  PointIndex pi = bel.PNum(k);
		  FuncDeriv_pts2d.Elem(pi) = new Point2d(ld.t1.functorMethod * (new mesh.Point(pi) - ld.sp1), ld.t2.functorMethod * (new mesh.Point(pi) - ld.sp1));
	  }
	  FuncDeriv_pts2d.Elem(gpi) = new Point2d(x(0), x(1));


	  vdir = new Vec2d(dir(0), dir(1));

	  hbad = bel.CalcJacobianBadnessDirDeriv(FuncDeriv_pts2d, lpi, vdir, ref hderiv);

	  deriv += hderiv;
	  badness += hbad;
	  }


	  return badness;
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual double Func(const Vector & x) const
	public override double Func(Vector x)
	{
	  Vector g = new Vector(x.Size());
	  return FuncGrad(x, ref g);
	}
  }





//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void MeshOptimize2d::GetNormalVector(int, const Point<3> & p, Vec<3> & nv) const

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void MeshOptimize2d::GetNormalVector(int surfind, const Point<3> & p, PointGeomInfo & gi, Vec<3> & n) const
}
