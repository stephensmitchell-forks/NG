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
#if SOLIDGEOM
#endif


namespace netgen
{


//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: double MinFunctionSum::Func(const Vector & x) const

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void MinFunctionSum::Grad(const Vector & x, Vector & g) const


//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: double MinFunctionSum::FuncGrad(const Vector & x, Vector & g) const

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: double MinFunctionSum::FuncDeriv(const Vector & x, const Vector & dir, double & deriv) const

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: double MinFunctionSum::GradStopping(const Vector & x) const



//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const MinFunction & MinFunctionSum::Function(int i) const



//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: double PointFunction1::Func(const Vector & vp) const


//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: double PointFunction1::FuncDeriv(const Vector & x, const Vector & dir, double & deriv) const


//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: double PointFunction1::FuncGrad(const Vector & x, Vector & g) const

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: double PointFunction1::GradStopping(const Vector & x) const


  /* Cheap Functional depending of inner point inside triangular surface */

  // is it used ????
  public class CheapPointFunction1 : MinFunction
  {
	private Mesh.T_POINTS points;
	private readonly Array<INDEX_3> faces;
	private DenseMatrix m = new DenseMatrix();
	private double h;
	public CheapPointFunction1(Mesh.T_POINTS apoints, Array<INDEX_3> afaces, double ah)
	{
		this.points = apoints;
		this.faces = new Array<INDEX_3>(afaces);
	  h = ah;


	  int nf = faces.Size();

	  m.SetSize(nf, 4);

	  for (int i = 1; i <= nf; i++)
	  {
	  Point3d p1 = points[new PointIndex(faces.Get(i).I1())];
	  Point3d p2 = points[new PointIndex(faces.Get(i).I2())];
	  Point3d p3 = points[new PointIndex(faces.Get(i).I3())];
	  Vec3d v1 = new Vec3d(p1, p2);
	  Vec3d v2 = new Vec3d(p1, p3);
	  Vec3d n = new Vec3d();
	  netgen.GlobalMembers.Cross(v1, v2, n);
	  n /= n.Length();

	  m.Elem(i, 1) = n.X();
	  m.Elem(i, 2) = n.Y();
	  m.Elem(i, 3) = n.Z();
	  m.Elem(i, 4) = - (n.X() * p1.X() + n.Y() * p1.Y() + n.Z() * p1.Z());
	  }
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual double Func(const Vector & vp) const
	public override double Func(Vector vp)
	{

	  /*
	    int j;
	    double badness = 0;
	    Point3d pp(vp.Get(1), vp.Get(2), vp.Get(3));
  
	    for (j = 1; j <= faces.Size(); j++)
	    {
	    const INDEX_3 & el = faces.Get(j);
  
	    double bad = CalcTetBadness (points.Get(el.I1()),
	    points.Get(el.I3()),
	    points.Get(el.I2()),
	    pp, 0);
	    badness += bad;
	    }
	  */

	  int i;
	  double badness = 0;
	  VectorMem < 4> hv;
	  Vector res = new Vector(m.Height());

	  for (i = 0;i < 3; i++)
	  {
		hv(i) = vp(i);
	  }
	  hv(3) = 1;
	  m.Mult(hv, res);

	  for (i = 1; i <= res.Size(); i++)
	  {
	  if (res(i - 1) < 1e-10)
	  {
		badness += 1e24;
	  }
	  else
	  {
		badness += 1 / res(i - 1);
	  }
	  }

	  return badness;
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual double FuncGrad(const Vector & x, Vector & g) const
	public override double FuncGrad(Vector x, Vector g)
	{
	  VectorMem < 3> hx;
	  double eps = 1e-6;

	  hx = x;
	  for (int i = 0; i < 3; i++)
	  {
	  hx(i) = x(i) + eps * h;
	  double fr = Func(hx);
	  hx(i) = x(i) - eps * h;
	  double fl = Func(hx);
	  hx(i) = x(i);

	  g(i) = (fr - fl) / (2 * eps * h);
	  }

	  return Func(x);
	}
  }














  /* ************* PointFunction **************************** */


  public class PointFunction : System.IDisposable
  {
	public Mesh.T_POINTS points;
	public readonly Array<Element, 0, uint> elements;
	public TABLE<int,PointIndex.BASE> elementsonpoint = new TABLE<int,PointIndex.BASE>();
	public readonly MeshingParameters mp;
	public PointIndex actpind = new PointIndex();
	public double h;

	public PointFunction(Mesh.T_POINTS apoints, Array<Element, 0, uint> aelements, MeshingParameters amp)
	{
		this.points = apoints;
		this.elements = new Array<Element, 0, uint>(aelements);
		this.elementsonpoint = new TABLE<int,PointIndex.BASE>(apoints.Size());
		this.mp = new netgen.MeshingParameters(amp);
	  for (int i = 0; i < elements.Size(); i++)
	  {
		if (elements[i].NP() == 4)
		{
		  for (int j = 0; j < elements[i].NP(); j++)
		  {
			elementsonpoint.Add(elements[i][j], i);
		  }
		}
	  }
	}

	public virtual void Dispose()
	{
		;
	}
	public virtual void SetPointIndex(PointIndex aactpind)
	{
	  actpind.CopyFrom(aactpind);
	}

	public void SetLocalH(double ah)
	{
		h = ah;
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: double GetLocalH() const
	public double GetLocalH()
	{
		return h;
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual double PointFunctionValue(const Point<3> & pp) const
	public virtual double PointFunctionValue(Point < 3> pp)
	{
	  double badness;
	  Point < 3> hp;

	  badness = 0;

	  hp = points[actpind];
	  points[actpind] = Point < 3> (pp.functorMethod);

	  for (int j = 0; j < elementsonpoint[actpind].Size(); j++)
	  {
		  Element el = elements[elementsonpoint[actpind][j]];
	  badness += CalcTetBadness(points[el[0]], points[el[1]], points[el[2]], points[el[3]], -1, mp);
	  }

	  points[actpind] = Point < 3> (hp);
	  return badness;
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual double PointFunctionValueGrad(const Point<3> & pp, Vec<3> & grad) const
	public virtual double PointFunctionValueGrad(Point < 3> pp, ref Vec < 3> grad)
	{
	  double f = 0;

	  Point < 3> hp = points[actpind];
	  Vec < 3> vgradi, vgrad(0,0,0);
	  points[actpind] = Point < 3> (pp.functorMethod);

	  for (int j = 0; j < elementsonpoint[actpind].Size(); j++)
	  {
		  Element el = elements[elementsonpoint[actpind][j]];
	  for (int k = 0; k < 4; k++)
	  {
		if (el[k] == actpind)
		{
			f += CalcTetBadnessGrad(points[el[0]], points[el[1]], points[el[2]], points[el[3]], -1, k + 1, vgradi, mp);

				vgrad += vgradi;
		}
	  }
	  }

	  points[actpind] = Point < 3> (hp);

	  grad = vgrad;
	  return f;
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual double PointFunctionValueDeriv(const Point<3> & pp, const Vec<3> & dir, double & deriv) const
	public virtual double PointFunctionValueDeriv(Point < 3> pp, Vec < 3> dir, ref double deriv)
	{
	  Vec < 3> vgradi, vgrad(0,0,0);

	  Point < 3> hp = points[actpind];
	  points[actpind] = pp.functorMethod;
	  double f = 0;

	  for (int j = 0; j < elementsonpoint[actpind].Size(); j++)
	  {
		  Element el = elements[elementsonpoint[actpind][j]];

	  for (int k = 1; k <= 4; k++)
	  {
		if (el.PNum(k) == actpind)
		{
			f += CalcTetBadnessGrad(points[el.PNum(1)], points[el.PNum(2)], points[el.PNum(3)], points[el.PNum(4)], -1, k, vgradi, mp);

			vgrad += vgradi;
		}
	  }
	  }

	  points[actpind] = Point < 3> (hp);
	  deriv = dir.functorMethod * vgrad;
	  return f;
	}

	public int MovePointToInner()
	{
	  // try point movement
	  Array<Element2d> faces = new Array<Element2d>();

	  for (int j = 0; j < elementsonpoint[actpind].Size(); j++)
	  {
	  Element el = elements[elementsonpoint[actpind][j]];

	  for (int k = 1; k <= 4; k++)
	  {
		if (el.PNum(k) == actpind)
		{
			Element2d face = new Element2d(ELEMENT_TYPE.TRIG);
			el.GetFace(k, face);
			netgen.GlobalMembers.Swap(ref face.PNum(2), ref face.PNum(3));
			faces.Append(face);
		}
	  }
	  }

	  Point3d hp = new Point3d();
	  int hi = netgen.GlobalMembers.FindInnerPoint(points, faces, hp);
	  if (hi != 0)
	  {
	  // cout << "inner point found" << endl;
	  points[actpind] = Point < 3> (hp);
	  }
	  else
	  {
		;
	  }
	  //      cout << "no inner point found" << endl;

	  /*
	  Point3d hp2;
	  int hi2 = FindInnerPoint (points, faces, hp2);
	  if (hi2)
	    {
	  cout << "new: inner point found" << endl;
	    }
	  else
	    cout << "new: no inner point found" << endl;
  
	  (*testout) << "hi(orig) = " << hi << ", hi(new) = " << hi2;
	  if (hi != hi2) (*testout) << "hi different" << endl;
	  */

	  return hi;
	}
  }






  public class CheapPointFunction : PointFunction
  {
	private DenseMatrix m = new DenseMatrix();
	public CheapPointFunction(Mesh.T_POINTS apoints, Array<Element, 0, uint> aelements, MeshingParameters amp) : base(apoints, aelements, amp)
	{
	  ;
	}

	public override void SetPointIndex(PointIndex aactpind)
	{
	  actpind.CopyFrom(aactpind);

	  int ne = elementsonpoint[actpind].Size();
	  int i;
	  int j;
	  PointIndex pi1 = new PointIndex();
	  PointIndex pi2 = new PointIndex();
	  PointIndex pi3 = new PointIndex();

	  m.SetSize(ne, 4);

	  for (i = 0; i < ne; i++)
	  {
	  pi1 = 0;
	  pi2 = 0;
	  pi3 = 0;

	  Element el = elements[elementsonpoint[actpind][i]];
	  for (j = 1; j <= 4; j++)
	  {
		if (el.PNum(j) != actpind)
		{
			pi3.CopyFrom(pi2);
			pi2.CopyFrom(pi1);
			pi1.CopyFrom(el.PNum(j));
		}
	  }

	  Point3d p1 = points[pi1];
	  Vec3d v1 = new Vec3d(p1, points[pi2]);
	  Vec3d v2 = new Vec3d(p1, points[pi3]);
	  Vec3d n = new Vec3d();
	  netgen.GlobalMembers.Cross(v1, v2, n);
	  n /= n.Length();

	  Vec3d v = new Vec3d(p1, points[actpind]);
	  double c = v * n;

	  if (c < 0)
	  {
		n *= -1;
	  }

	  // n is inner normal

	  m.Elem(i + 1, 1) = n.X();
	  m.Elem(i + 1, 2) = n.Y();
	  m.Elem(i + 1, 3) = n.Z();
	  m.Elem(i + 1, 4) = - (n.X() * p1.X() + n.Y() * p1.Y() + n.Z() * p1.Z());
	  }
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual double PointFunctionValue(const Point<3> & pp) const
	public override double PointFunctionValue(Point < 3> pp)
	{
	  VectorMem < 4> p4;
	  Vector di = new Vector();
	  int n = m.Height();

	  p4(0) = pp.functorMethod(0);
	  p4(1) = pp.functorMethod(1);
	  p4(2) = pp.functorMethod(2);
	  p4(3) = 1;

	  di.SetSize(n);
	  m.Mult(p4, di);

	  double sum = 0;
	  for (int i = 0; i < n; i++)
	  {
	  if (di(i) > 0)
	  {
		sum += 1 / di(i);
	  }
	  else
	  {
		return 1e16;
	  }
	  }
	  return sum;
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual double PointFunctionValueGrad(const Point<3> & pp, Vec<3> & grad) const
	public override double PointFunctionValueGrad(Point < 3> pp, ref Vec < 3> grad)
	{
	  VectorMem < 4> p4;
	  Vector di = new Vector();

	  int n = m.Height();

	  p4(0) = pp.functorMethod(0);
	  p4(1) = pp.functorMethod(1);
	  p4(2) = pp.functorMethod(2);
	  p4(3) = 1;

	  di.SetSize(n);
	  m.Mult(p4, di);

	  double sum = 0;
	  grad = 0;
	  for (int i = 0; i < n; i++)
	  {
	  if (di(i) > 0)
	  {
		  double idi = 1 / di(i);
		  sum += idi;
		  grad.functorMethod(0) -= idi * idi * m.functorMethod(i, 0);
		  grad.functorMethod(1) -= idi * idi * m.functorMethod(i, 1);
		  grad.functorMethod(2) -= idi * idi * m.functorMethod(i, 2);
	  }
	  else
	  {
		  return 1e16;
	  }
	  }
	  return sum;
	}
  }








  public class Opti3FreeMinFunction : MinFunction
  {
	private readonly PointFunction pf;
	private Point < 3> sp1;

	public Opti3FreeMinFunction(PointFunction apf)
	{
		this.pf = new netgen.PointFunction(apf);
	  ;
	}

	public void SetPoint(Point < 3> asp1)
	{
		sp1 = asp1.functorMethod;
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual double Func(const Vector & x) const
	public override double Func(Vector x)
	{
	  Point < 3> pp;
	  for (int j = 0; j < 3; j++)
	  {
		pp(j) = sp1(j) + x(j);
	  }
	  return pf.PointFunctionValue(pp);
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual double FuncGrad(const Vector & x, Vector & grad) const
	public override double FuncGrad(Vector x, Vector grad)
	{
	  Vec < 3> vgrad;
	  Point < 3> pp;

	  for (int j = 0; j < 3; j++)
	  {
		pp(j) = sp1(j) + x(j);
	  }

	  double val = pf.PointFunctionValueGrad(pp, ref vgrad);

	  for (int j = 0; j < 3; j++)
	  {
		grad(j) = vgrad(j);
	  }

	  return val;
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual double FuncDeriv(const Vector & x, const Vector & dir, double & deriv) const
	public override double FuncDeriv(Vector x, Vector dir, ref double deriv)
	{
	  Point < 3> pp;

	  for (int j = 0; j < 3; j++)
	  {
		pp(j) = sp1(j) + x(j);
	  }

	  Vec < 3> vdir;
	  for (int j = 0; j < 3; j++)
	  {
		vdir(j) = dir(j);
	  }

	  return pf.PointFunctionValueDeriv(pp, vdir, ref deriv);
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual double GradStopping(const Vector & x) const
	public override double GradStopping(Vector x)
	{
	  double f = Func(x);
	  return 1e-3 * f / pf.GetLocalH();
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual void ApproximateHesse(const Vector & x, DenseMatrix & hesse) const
	public override void ApproximateHesse(Vector x, DenseMatrix hesse)
	{
	  int n = x.Size();

	  Vector hx = new Vector();
	  hx.SetSize(n);

	  double eps = 1e-8;
	  double f; //, f12, f21
	  double f11;
	  double f22;

	  f = Func(x);

	  for (int i = 1; i <= n; i++)
	  {
	  for (int j = 1; j < i; j++)
	  {
		  /*
		    hx = x;
		    hx.Elem(i) = x.Get(i) + eps;
		    hx.Elem(j) = x.Get(j) + eps;
		    f11 = Func(hx);
		    hx.Elem(i) = x.Get(i) + eps;
		    hx.Elem(j) = x.Get(j) - eps;
		    f12 = Func(hx);
		    hx.Elem(i) = x.Get(i) - eps;
		    hx.Elem(j) = x.Get(j) + eps;
		    f21 = Func(hx);
		    hx.Elem(i) = x.Get(i) - eps;
		    hx.Elem(j) = x.Get(j) - eps;
		    f22 = Func(hx);
		  */
		  hesse.Elem(i, j) = hesse.Elem(j, i) = 0;
		  //	    (f11 + f22 - f12 - f21) / (2 * eps * eps);
	  }

	  hx.CopyFrom(x);
	  hx(i - 1) = x(i - 1) + eps;
	  f11 = Func(hx);
	  hx(i - 1) = x(i - 1) - eps;
	  f22 = Func(hx);

	  hesse.Elem(i, i) = (f11 + f22 - 2 * f) / (eps * eps) + 1e-12;
	  }
	}
  }






#if SOLIDGEOM
  public class Opti3SurfaceMinFunction : MinFunction
  {
	private readonly PointFunction pf;
	private Point3d sp1 = new Point3d();
	private readonly Surface surf;
	private Vec3d t1 = new Vec3d();
	private Vec3d t2 = new Vec3d();

	public Opti3SurfaceMinFunction(PointFunction apf) : base()
	{
		this.pf = new netgen.PointFunction(apf);
	  ;
	}

	public void SetPoint(Surface asurf, Point3d asp1)
	{
	  Vec3d n = new Vec3d();
	  sp1.CopyFrom(asp1);
	  surf = asurf;

	  Vec < 3> hn;
	  surf.GetNormalVector(sp1, hn);
	  n = hn;

	  n.GetNormal(t1);
	  t1 /= t1.Length();
	  t2 = netgen.GlobalMembers.Cross(n, t1);
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void CalcNewPoint(const Vector & x, Point3d & np) const
	public void CalcNewPoint(Vector x, ref Point3d np)
	{
	  np.X() = sp1.X() + x.Get(1) * t1.X() + x.Get(2) * t2.X();
	  np.Y() = sp1.Y() + x.Get(1) * t1.Y() + x.Get(2) * t2.Y();
	  np.Z() = sp1.Z() + x.Get(1) * t1.Z() + x.Get(2) * t2.Z();

	  Point < 3> hnp = np;
	  surf.Project(hnp);
	  np = hnp;
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual double Func(const Vector & x) const
	public override double Func(Vector x)
	{
	  Point3d pp1 = new Point3d();

	  CalcNewPoint(x, ref pp1);
	  return pf.PointFunctionValue(pp1);
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual double FuncGrad(const Vector & x, Vector & grad) const
	public override double FuncGrad(Vector x, Vector grad)
	{
	  Vec3d n = new Vec3d();
	  Vec3d vgrad = new Vec3d();
	  Point3d pp1 = new Point3d();
	  VectorMem < 3> freegrad;

	  CalcNewPoint(x, ref pp1);

	  double badness = pf.PointFunctionValueGrad(pp1, ref freegrad);
	  vgrad.X() = freegrad.Get(1);
	  vgrad.Y() = freegrad.Get(2);
	  vgrad.Z() = freegrad.Get(3);

	  Vec < 3> hn;
	  surf.GetNormalVector(pp1, hn);
	  n = hn;

	  vgrad -= (vgrad * n) * n;

	  grad.Elem(1) = vgrad * t1;
	  grad.Elem(2) = vgrad * t2;

	  return badness;
	}
  }
#endif








#if SOLIDGEOM
  public class Opti3EdgeMinFunction : MinFunction
  {
	private readonly PointFunction pf;
	private Point3d sp1 = new Point3d();
	private readonly Surface surf1;
	private readonly Surface surf2;
	private Vec3d t1 = new Vec3d();

	public Opti3EdgeMinFunction(PointFunction apf) : base()
	{
		this.pf = new netgen.PointFunction(apf);
	  ;
	}

	public void SetPoint(Surface asurf1, Surface asurf2, Point3d asp1)
	{
	  Vec3d n1 = new Vec3d();
	  Vec3d n2 = new Vec3d();
	  sp1.CopyFrom(asp1);
	  surf1 = asurf1;
	  surf2 = asurf2;

	  Vec < 3> hn1, hn2;
	  surf1.GetNormalVector(sp1, hn1);
	  surf2.GetNormalVector(sp1, hn2);
	  n1 = hn1;
	  n2 = hn2;
	  t1 = netgen.GlobalMembers.Cross(n1, n2);
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void CalcNewPoint(const Vector & x, Point3d & np) const
	public void CalcNewPoint(Vector x, ref Point3d np)
	{
	np.X() = sp1.X() + x.Get(1) * t1.X();
	np.Y() = sp1.Y() + x.Get(1) * t1.Y();
	np.Z() = sp1.Z() + x.Get(1) * t1.Z();
	Point < 3> hnp = np;
	ProjectToEdge(surf1, surf2, hnp);
	np = hnp;
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual double FuncGrad(const Vector & x, Vector & grad) const
	public override double FuncGrad(Vector x, Vector grad)
	{
	  Vec3d n1 = new Vec3d();
	  Vec3d n2 = new Vec3d();
	  Vec3d v1 = new Vec3d();
	  Vec3d vgrad = new Vec3d();
	  Point3d pp1 = new Point3d();
	  double badness;
	  VectorMem < 3> freegrad;

	  CalcNewPoint(x, ref pp1);


	  badness = pf.PointFunctionValueGrad(pp1, ref freegrad);

	  vgrad.X() = freegrad.Get(1);
	  vgrad.Y() = freegrad.Get(2);
	  vgrad.Z() = freegrad.Get(3);

	  Vec < 3> hn1, hn2;
	  surf1.GetNormalVector(pp1, hn1);
	  surf2.GetNormalVector(pp1, hn2);
	  n1 = hn1;
	  n2 = hn2;

	  v1 = netgen.GlobalMembers.Cross(n1, n2);
	  v1 /= v1.Length();

	  grad.Elem(1) = (vgrad * v1) * (t1 * v1);
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
#endif











/* ************* JacobianPointFunction **************************** */




// class JacobianPointFunction : public MinFunction
// {
// public:
//   Mesh::T_POINTS & points;
//   const Mesh::T_VOLELEMENTS & elements;
//   TABLE<INDEX> elementsonpoint;
//   PointIndex actpind;

// public:
//   JacobianPointFunction (Mesh::T_POINTS & apoints, 
// 			 const Mesh::T_VOLELEMENTS & aelements);

//   virtual void SetPointIndex (PointIndex aactpind);
//   virtual double Func (const Vector & x) const;
//   virtual double FuncGrad (const Vector & x, Vector & g) const;
//   virtual double FuncDeriv (const Vector & x, const Vector & dir, double & deriv) const;
// };





//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: double JacobianPointFunction::Func(const Vector & v) const





//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: double JacobianPointFunction::FuncGrad(const Vector & x, Vector & g) const


//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: double JacobianPointFunction::FuncDeriv(const Vector & x, const Vector & dir, double & deriv) const










#if SOLIDGEOMxxxx
#endif




//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
//C++ TO C# CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
Timer ImproveMesh_t("Mesh::ImproveMesh");





// Improve Condition number of Jacobian, any elements  




// Improve Condition number of Jacobian, any elements  




}
