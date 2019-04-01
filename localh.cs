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
/* File:   localh.hh                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   29. Jan. 97                                                    */
/**************************************************************************/


namespace netgen
{


  /// box for grading
  public class GradingBox
  {
	/// xmid
	private float[] xmid = new float[3];
	/// half edgelength
	private float h2;
	///
	private GradingBox[] childs = Arrays.InitializeWithDefaultInstances<GradingBox>(8);
	///
	private GradingBox father;
	///
	private double hopt;
	///

//C++ TO C# CONVERTER NOTE: Classes must be named in C#, so the following class has been named by the converter:
	public class AnonymousClass
	{
//C++ TO C# CONVERTER TODO TASK: C# does not allow bit fields:
	  public uint cutboundary:1;
//C++ TO C# CONVERTER TODO TASK: C# does not allow bit fields:
	  public uint isinner:1;
//C++ TO C# CONVERTER TODO TASK: C# does not allow bit fields:
	  public uint oldcell:1;
//C++ TO C# CONVERTER TODO TASK: C# does not allow bit fields:
	  public uint pinner:1;
	}
	public AnonymousClass flags = new AnonymousClass();

	///
	public GradingBox(double[] ax1, double[] ax2)
	{
	  h2 = 0.5 * (ax2[0] - ax1[0]);
	  for (int i = 0; i < 3; i++)
	  {
		xmid[i] = 0.5 * (ax1[i] + ax2[i]);
	  }

	  for (int i = 0; i < 8; i++)
	  {
		childs[i] = null;
	  }
	  father = null;

	  flags.cutboundary = 0;
	  flags.isinner = 0;
	  flags.oldcell = 0;
	  flags.pinner = 0;

	  hopt = 2 * h2;
	}

	///
	public void DeleteChilds()
	{
	  for (int i = 0; i < 8; i++)
	  {
		if (childs[i] != null)
		{
		childs[i].DeleteChilds();
		childs[i] = null;
		childs[i] = null;
		}
	  }
	}

	///

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: Point<3> PMid() const
	public Point < 3> PMid()
	{
		return Point < 3> (xmid[0], xmid[1], xmid[2]);
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: double H2() const
	public double H2()
	{
		return h2;
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool HasChilds() const
	public bool HasChilds()
	{
	  for (int i = 0; i < 8; i++)
	  {
		if (childs[i] != null)
		{
			return true;
		}
	  }
	  return false;
	}

//C++ TO C# CONVERTER TODO TASK: C# has no concept of a 'friend' class:
//	friend class LocalH;

	public static BlockAllocator ball = new BlockAllocator(sizeof(GradingBox));
//C++ TO C# CONVERTER TODO TASK: The new operator cannot be overloaded in C#:
	public static object operator new(uint UnnamedParameter)
	{
	  return ball.Alloc();
	}

//C++ TO C# CONVERTER TODO TASK: The delete operator cannot be overloaded in C#:
	public static void operator delete(object p)
	{
	  ball.Free(p);
	}
  }




  /**
     Control of 3D mesh grading
  */
  public class LocalH : System.IDisposable
  {
	///
	private GradingBox root;
	///
	private double grading;
	///
	private Array<GradingBox> boxes = new Array<GradingBox>();
	///
	private Box < 3> boundingbox;
	/// octree or quadtree
	private int dimension;
	///
	public LocalH(Point < 3> pmin, Point < 3> pmax, double agrading, int adimension = 3)
	{
		this.dimension = adimension;
	  double[] x1 = new double[3];
	  double[] x2 = new double[3];
	  double hmax;

	  boundingbox = Box < 3> (pmin.functorMethod, pmax.functorMethod);
	  grading = agrading;

	  // a small enlargement, non-regular points
	  double val = 0.0879;
	  for (int i = 0; i < dimension; i++)
	  {
	  x1[i] = (1 + val * (i + 1)) * pmin.functorMethod(i) - val * (i + 1) * pmax.functorMethod(i);
	  x2[i] = 1.1 * pmax.functorMethod(i) - 0.1 * pmin.functorMethod(i);
	  }
	  for (int i = dimension; i < 3; i++)
	  {
		x1[i] = x2[i] = 0;
	  }

	  hmax = x2[0] - x1[0];
	  for (int i = 1; i < dimension; i++)
	  {
		hmax = netgen.GlobalMembers.max2(x2[i] - x1[i], hmax);
	  }

	  for (int i = 0; i < dimension; i++)
	  {
		x2[i] = x1[i] + hmax;
	  }

	  root = new GradingBox(x1, x2);
	  boxes.Append(root);
	}

	///
	public LocalH(Box < 3> box, double grading, int adimension = 3) : this(box.PMin(), box.PMax(), grading, adimension)
	{
		;
	}
	///
	public void Dispose()
	{
	  root.DeleteChilds();
	  if (root != null)
	  {
		  root.Dispose();
	  }
	}

	///
	public void Delete()
	{
	  root.DeleteChilds();
	}

	///
	public void SetGrading(double agrading)
	{
		grading = agrading;
	}
	///
	public void SetH(Point < 3> p, double h)
	{
	  if (dimension == 2)
	  {
		  if (ngsimd.GlobalMembers.fabs(p.functorMethod(0) - root.xmid[0]) > root.h2 || ngsimd.GlobalMembers.fabs(p.functorMethod(1) - root.xmid[1]) > root.h2)
		  {
			return;
		  }

		  if (GetH(p.functorMethod) <= 1.2 * h)
		  {
			  return;
		  }

		  GradingBox box = root;
		  GradingBox nbox = root;
		  GradingBox ngb;
		  int childnr;
		  double[] x1 = new double[3];
		  double[] x2 = new double[3];

		  while (nbox != null)
		  {
			  box = nbox;
			  childnr = 0;
			  if (p.functorMethod(0) > box.xmid[0])
			  {
				  childnr += 1;
			  }
			  if (p.functorMethod(1) > box.xmid[1])
			  {
				  childnr += 2;
			  }
			  nbox = box.childs[childnr];
		  };

		  while (2 * box.h2 > h)
		  {
			  childnr = 0;
			  if (p.functorMethod(0) > box.xmid[0])
			  {
				  childnr += 1;
			  }
			  if (p.functorMethod(1) > box.xmid[1])
			  {
				  childnr += 2;
			  }

			  double h2 = box.h2;
			  if ((childnr & 1) != 0)
			  {
				  x1[0] = box.xmid[0];
				  x2[0] = x1[0] + h2; // box->x2[0];
			  }
			  else
			  {
				  x2[0] = box.xmid[0];
				  x1[0] = x2[0] - h2; // box->x1[0];
			  }

			  if ((childnr & 2) != 0)
			  {
				  x1[1] = box.xmid[1];
				  x2[1] = x1[1] + h2; // box->x2[1];
			  }
			  else
			  {
				  x2[1] = box.xmid[1];
				  x1[1] = x2[1] - h2; // box->x1[1];
			  }
			  x1[2] = x2[2] = 0;

			  ngb = new GradingBox(x1, x2);
			  box.childs[childnr] = ngb;
			  ngb.father = box;

			  boxes.Append(ngb);
			  box = box.childs[childnr];
		  }

		  box.hopt = h;

		  double hbox = 2 * box.h2; // box->x2[0] - box->x1[0];
		  double hnp = h + grading * hbox;

		  Point < 3> np;
		  for (int i = 0; i < 2; i++)
		  {
			  np = p.functorMethod;
			  np(i) = p.functorMethod(i) + hbox;
			  SetH(np, hnp);

			  np(i) = p.functorMethod(i) - hbox;
			  SetH(np, hnp);
		  }

	  }

	  else
	  {
		  if (ngsimd.GlobalMembers.fabs(p.functorMethod(0) - root.xmid[0]) > root.h2 || ngsimd.GlobalMembers.fabs(p.functorMethod(1) - root.xmid[1]) > root.h2 || ngsimd.GlobalMembers.fabs(p.functorMethod(2) - root.xmid[2]) > root.h2)
		  {
			return;
		  }

		  if (GetH(p.functorMethod) <= 1.2 * h)
		  {
			  return;
		  }

		  GradingBox box = root;
		  GradingBox nbox = root;
		  GradingBox ngb;
		  int childnr;
		  double[] x1 = new double[3];
		  double[] x2 = new double[3];

		  while (nbox != null)
		  {
			  box = nbox;
			  childnr = 0;
			  if (p.functorMethod(0) > box.xmid[0])
			  {
				  childnr += 1;
			  }
			  if (p.functorMethod(1) > box.xmid[1])
			  {
				  childnr += 2;
			  }
			  if (p.functorMethod(2) > box.xmid[2])
			  {
				  childnr += 4;
			  }
			  nbox = box.childs[childnr];
		  };


		  while (2 * box.h2 > h)
		  {
			  childnr = 0;
			  if (p.functorMethod(0) > box.xmid[0])
			  {
				  childnr += 1;
			  }
			  if (p.functorMethod(1) > box.xmid[1])
			  {
				  childnr += 2;
			  }
			  if (p.functorMethod(2) > box.xmid[2])
			  {
				  childnr += 4;
			  }

			  double h2 = box.h2;
			  if ((childnr & 1) != 0)
			  {
				  x1[0] = box.xmid[0];
				  x2[0] = x1[0] + h2; // box->x2[0];
			  }
			  else
			  {
				  x2[0] = box.xmid[0];
				  x1[0] = x2[0] - h2; // box->x1[0];
			  }

			  if ((childnr & 2) != 0)
			  {
				  x1[1] = box.xmid[1];
				  x2[1] = x1[1] + h2; // box->x2[1];
			  }
			  else
			  {
				  x2[1] = box.xmid[1];
				  x1[1] = x2[1] - h2; // box->x1[1];
			  }

			  if ((childnr & 4) != 0)
			  {
				  x1[2] = box.xmid[2];
				  x2[2] = x1[2] + h2; // box->x2[2];
			  }
			  else
			  {
				  x2[2] = box.xmid[2];
				  x1[2] = x2[2] - h2; // box->x1[2];
			  }

			  ngb = new GradingBox(x1, x2);
			  box.childs[childnr] = ngb;
			  ngb.father = box;

			  boxes.Append(ngb);
			  box = box.childs[childnr];
		  }

		  box.hopt = h;

		  double hbox = 2 * box.h2; // box->x2[0] - box->x1[0];
		  double hnp = h + grading * hbox;

		  Point < 3> np;
		  for (int i = 0; i < 3; i++)
		  {
			  np = p.functorMethod;
			  np(i) = p.functorMethod(i) + hbox;
			  SetH(np, hnp);

			  np(i) = p.functorMethod(i) - hbox;
			  SetH(np, hnp);
		  }
	  }
	}

	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: double GetH(Point<3> x) const
	public double GetH(Point < 3> x)
	{
	  GradingBox box = root;
	  if (dimension == 2)
	  {
		  while (true)
		  {
			  int childnr = 0;
			  if (x.functorMethod(0) > box.xmid[0])
			  {
				  childnr += 1;
			  }
			  if (x.functorMethod(1) > box.xmid[1])
			  {
				  childnr += 2;
			  }

			  if (box.childs[childnr])
			  {
				box = box.childs[childnr];
			  }
			  else
			  {
				return box.hopt;
			  }
		  }
	  }
	  else
	  {
		  while (true)
		  {
			  int childnr = 0;
			  if (x.functorMethod(0) > box.xmid[0])
			  {
				  childnr += 1;
			  }
			  if (x.functorMethod(1) > box.xmid[1])
			  {
				  childnr += 2;
			  }
			  if (x.functorMethod(2) > box.xmid[2])
			  {
				  childnr += 4;
			  }

			  if (box.childs[childnr])
			  {
				box = box.childs[childnr];
			  }
			  else
			  {
				return box.hopt;
			  }
		  }
	  }

	}

	/// minimal h in box (pmin, pmax)

	/// minimal h in box (pmin, pmax)
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: double GetMinH(Point<3> pmin, Point<3> pmax) const
	public double GetMinH(Point < 3> pmin, Point < 3> pmax)
	{
	  Point < 3> pmin2, pmax2;
	  for (int j = 0; j < 3; j++)
	  {
		if (pmin.functorMethod(j) < pmax.functorMethod(j))
		{
		  pmin2(j) = pmin.functorMethod(j);
		  pmax2(j) = pmax.functorMethod(j);
		}
		else
		{
		  pmin2(j) = pmax.functorMethod(j);
		  pmax2(j) = pmin.functorMethod(j);
		}
	  }

	  return GetMinHRec(pmin2, pmax2, root);
	}

	/// mark boxes intersecting with boundary-box
	// void CutBoundary (const Point3d & pmin, const Point3d & pmax)
	// { CutBoundaryRec (pmin, pmax, root); }

	public void CutBoundary(Box < 3> box)
	{
		CutBoundaryRec(box.PMin(), box.PMax(), root);
	}

	/// find inner boxes
	private delegate int testinnerDelegate(Point3d p1);

//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//	void FindInnerBoxes(AdFront3 adfront, testinnerDelegate testinner);

	private delegate int testinnerDelegate(Point < 2> p1);

//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//	void FindInnerBoxes(AdFront2 adfront, testinnerDelegate testinner);


	/// clears all flags 
	public void ClearFlags()
	{
		ClearFlagsRec(root);
	}

	/// widen refinement zone
	public void WidenRefinement()
	{
	  for (int i = 0; i < boxes.Size(); i++)
	  {
	  double h = boxes[i].hopt;
	  Point3d c = boxes[i].PMid();

	  for (int i1 = -1; i1 <= 1; i1++)
	  {
		for (int i2 = -1; i2 <= 1; i2++)
		{
		  for (int i3 = -1; i3 <= 1; i3++)
		  {
			SetH(new Point3d(c.X() + i1 * h, c.Y() + i2 * h, c.Z() + i3 * h), 1.001 * h);
		  }
		}
	  }
	  }
	}

	/// get points in inner elements
	public void GetInnerPoints(Array<Point < 3>> points)
	{
	  if (dimension == 2)
	  {
		  for (int i = 0; i < boxes.Size(); i++)
		  {
			if (boxes[i].flags.isinner && boxes[i].HasChilds())
			{
			  points.Append(boxes[i].PMid());
			}
		  }
	  }
	  else
	  {
		  for (int i = 0; i < boxes.Size(); i++)
		  {
			if (boxes[i].flags.isinner)
			{
			  points.Append(boxes[i].PMid());
			}
		  }
	  }

	}

	/// get points in outer closure
	public void GetOuterPoints(Array<Point < 3>> points)
	{
	  for (int i = 0; i < boxes.Size(); i++)
	  {
		if (!boxes[i].flags.isinner && !boxes[i].flags.cutboundary)
		{
	  points.Append(boxes[i].PMid());
		}
	  }
	}

	///
	public void Convexify()
	{
	  ConvexifyRec(root);
	}

	///
	public int GetNBoxes()
	{
		return boxes.Size();
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const Box<3> & GetBoundingBox() const
	public Box < 3> GetBoundingBox()
	{
		return boundingbox;
	}
	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void PrintMemInfo(ostream & ost) const
	public void PrintMemInfo(ostream ost)
	{
	  ost << "LocalH: " << boxes.Size() << " boxes of " << sizeof(GradingBox) << " bytes = " << boxes.Size() * sizeof(GradingBox) << " bytes" << "\n";
	}

	/// 
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: double GetMinHRec(const Point3d & pmin, const Point3d & pmax, const GradingBox * box) const
	private double GetMinHRec(Point3d pmin, Point3d pmax, GradingBox box)
	{
	  if (dimension == 2)
	  {
		  double h2 = box.h2;
		  if (pmax.X() < box.xmid[0] - h2 || pmin.X() > box.xmid[0] + h2 || pmax.Y() < box.xmid[1] - h2 || pmin.Y() > box.xmid[1] + h2)
		  {
			return 1e8;
		  }

		  double hmin = 2 * box.h2; // box->x2[0] - box->x1[0];

		  for (int i = 0; i < 8; i++)
		  {
			if (box.childs[i])
			{
			  hmin = netgen.GlobalMembers.min2(hmin, GetMinHRec(pmin, pmax, box.childs[i]));
			}
		  }

		  return hmin;
	  }
	  else
	  {
		  double h2 = box.h2;
		  if (pmax.X() < box.xmid[0] - h2 || pmin.X() > box.xmid[0] + h2 || pmax.Y() < box.xmid[1] - h2 || pmin.Y() > box.xmid[1] + h2 || pmax.Z() < box.xmid[2] - h2 || pmin.Z() > box.xmid[2] + h2)
		  {
			return 1e8;
		  }

		  double hmin = 2 * box.h2; // box->x2[0] - box->x1[0];

		  for (int i = 0; i < 8; i++)
		  {
			if (box.childs[i])
			{
			  hmin = netgen.GlobalMembers.min2(hmin, GetMinHRec(pmin, pmax, box.childs[i]));
			}
		  }

		  return hmin;
	  }
	}

	///
	private void CutBoundaryRec(Point3d pmin, Point3d pmax, GradingBox box)
	{
	  double h2 = box.h2;
	  if (dimension == 2)
	  {
		  if (pmax.X() < box.xmid[0] - h2 || pmin.X() > box.xmid[0] + h2 || pmax.Y() < box.xmid[1] - h2 || pmin.Y() > box.xmid[1] + h2)
		  {
			return;
		  }
	  }
	  else
	  {
		  if (pmax.X() < box.xmid[0] - h2 || pmin.X() > box.xmid[0] + h2 || pmax.Y() < box.xmid[1] - h2 || pmin.Y() > box.xmid[1] + h2 || pmax.Z() < box.xmid[2] - h2 || pmin.Z() > box.xmid[2] + h2)
		  {
			return;
		  }
	  }

	  box.flags.cutboundary = 1;
	  for (int i = 0; i < 8; i++)
	  {
		if (box.childs[i])
		{
	  CutBoundaryRec(pmin, pmax, box.childs[i]);
		}
	  }
	}

	///
	private delegate int innerDelegate(Point3d p);

	private void FindInnerBoxesRec(innerDelegate inner, GradingBox box)
	{
	  if (box.flags.cutboundary)
	  {
	  for (int i = 0; i < 8; i++)
	  {
		if (box.childs[i])
		{
		  FindInnerBoxesRec(inner, box.childs[i]);
		}
	  }
	  }
	  else
	  {
	  if (inner(box.PMid()) != 0)
	  {
		SetInnerBoxesRec(box);
	  }
	  }
	}

	///
	private void FindInnerBoxesRec2(GradingBox box, AdFront3 adfront, Array<Box3d> faceboxes, Array<int> faceinds, int nfinbox)
	{
	  if (box == null)
	  {
		  return;
	  }

	  GradingBox father = box.father;

	  Point3d c = new Point3d(box.xmid[0], box.xmid[1], box.xmid[2]);
	  Vec3d v = new Vec3d(box.h2, box.h2, box.h2);
	  Box3d boxc = new Box3d(c - v, c + v);

	  Point3d fc = new Point3d(father.xmid[0], father.xmid[1], father.xmid[2]);
	  Vec3d fv = new Vec3d(father.h2, father.h2, father.h2);
	  Box3d fboxc = new Box3d(fc - fv, fc + fv);

	  Box3d boxcfc = new Box3d(c, fc);

	  ArrayMem<int, 100> faceused = new ArrayMem<int, 100>();
	  ArrayMem<int, 100> faceused2 = new ArrayMem<int, 100>();
	  ArrayMem<int, 100> facenotused = new ArrayMem<int, 100>();

	  /*
	  faceused.SetSize(0);
	  facenotused.SetSize(0);
	  faceused2.SetSize(0);
	  */

	  for (int j = 1; j <= nfinbox; j++)
	  {
	  //      adfront->GetFaceBoundingBox (faceinds.Get(j), facebox);
	  Box3d facebox = faceboxes.Get(faceinds.Get(j));

	  if (boxc.Intersect(facebox) != 0)
	  {
		faceused.Append(faceinds.Get(j));
	  }
	  else
	  {
		facenotused.Append(faceinds.Get(j));
	  }

	  if (boxcfc.Intersect(facebox) != 0)
	  {
		faceused2.Append(faceinds.Get(j));
	  }
	  }

	  for (int j = 1; j <= faceused.Size(); j++)
	  {
		faceinds.Elem(j) = faceused.Get(j);
	  }
	  for (int j = 1; j <= facenotused.Size(); j++)
	  {
		faceinds.Elem(j + faceused.Size()) = facenotused.Get(j);
	  }


	  if (!father.flags.cutboundary)
	  {
	  box.flags.isinner = father.flags.isinner;
	  box.flags.pinner = father.flags.pinner;
	  }
	  else
	  {
	  Point3d cf = new Point3d(father.xmid[0], father.xmid[1], father.xmid[2]);

	  if (father.flags.isinner)
	  {
		box.flags.pinner = 1;
	  }
	  else
	  {
		  if (adfront.SameSide(c, cf, faceused2) != 0)
		  {
			box.flags.pinner = father.flags.pinner;
		  }
		  else
		  {
			box.flags.pinner = 1 - father.flags.pinner;
		  }
	  }

	  if (box.flags.cutboundary)
	  {
		box.flags.isinner = 0;
	  }
	  else
	  {
		box.flags.isinner = box.flags.pinner;
	  }
	  }

	  // cout << "faceused: " << faceused.Size() << ", " << faceused2.Size() << ", " << facenotused.Size() << endl;

	  int nf = faceused.Size();
	  for (int i = 0; i < 8; i++)
	  {
		FindInnerBoxesRec2(box.childs[i], adfront, faceboxes, faceinds, nf);
	  }
	}



	private delegate int innerDelegate(Point < 2> p);

	private void FindInnerBoxesRec(innerDelegate inner, GradingBox box)
	{
	  if (box.flags.cutboundary)
	  {
	  for (int i = 0; i < 8; i++)
	  {
		if (box.childs[i])
		{
		  FindInnerBoxesRec(inner, box.childs[i]);
		}
	  }
	  }
	  else
	  {
	  Point < 2> p2d(box.PMid()(0), box.PMid()(1));
	  if (inner(p2d) != 0)
	  {
		SetInnerBoxesRec(box);
	  }
	  }
	}

	///
	private void FindInnerBoxesRec2(GradingBox box, AdFront2 adfront, Array<Box < 3>> faceboxes, Array<int> faceinds, int nfinbox)
	{
	  if (box == null)
	  {
		  return;
	  }

	  GradingBox father = box.father;

	  Point3d c = new Point3d(box.xmid[0], box.xmid[1], 0); // box->xmid[2]);
	  Vec3d v = new Vec3d(box.h2, box.h2, box.h2);
	  Box3d boxc = new Box3d(c - v, c + v);

	  Point3d fc = new Point3d(father.xmid[0], father.xmid[1], 0); // father->xmid[2]);
	  Vec3d fv = new Vec3d(father.h2, father.h2, father.h2);
	  Box3d fboxc = new Box3d(fc - fv, fc + fv);
	  Box3d boxcfc = new Box3d(c, fc);

	  ArrayMem<int, 100> faceused = new ArrayMem<int, 100>();
	  ArrayMem<int, 100> faceused2 = new ArrayMem<int, 100>();
	  ArrayMem<int, 100> facenotused = new ArrayMem<int, 100>();

	  for (int j = 0; j < nfinbox; j++)
	  {
	  //      adfront->GetFaceBoundingBox (faceinds.Get(j), facebox);
	  Box3d facebox = faceboxes[faceinds[j]];

	  if (boxc.Intersect(facebox) != 0)
	  {
		faceused.Append(faceinds[j]);
	  }
	  else
	  {
		facenotused.Append(faceinds[j]);
	  }

	  if (boxcfc.Intersect(facebox) != 0)
	  {
		faceused2.Append(faceinds[j]);
	  }
	  }

	  for (int j = 0; j < faceused.Size(); j++)
	  {
		faceinds[j] = faceused[j];
	  }
	  for (int j = 0; j < facenotused.Size(); j++)
	  {
		faceinds[j + faceused.Size()] = facenotused[j];
	  }

	  if (!father.flags.cutboundary)
	  {
	  box.flags.isinner = father.flags.isinner;
	  box.flags.pinner = father.flags.pinner;
	  }
	  else
	  {
	  Point3d cf = new Point3d(father.xmid[0], father.xmid[1], father.xmid[2]);

	  if (father.flags.isinner)
	  {
			  box.flags.pinner = 1;
	  }
	  else
	  {
		  Point < 2> c2d(c.X(), c.Y());
		  Point < 2> cf2d(cf.X(), cf.Y());
			  bool sameside = adfront.SameSide(c2d, cf2d, faceused2);
			  if (sameside)
			  {
			box.flags.pinner = father.flags.pinner;
			  }
		  else
		  {
			box.flags.pinner = 1 - father.flags.pinner;
		  }
	  }

	  if (box.flags.cutboundary)
	  {
		box.flags.isinner = 0;
	  }
	  else
	  {
		box.flags.isinner = box.flags.pinner;
	  }
	  }

	  // cout << "faceused: " << faceused.Size() << ", " << faceused2.Size() << ", " << facenotused.Size() << endl;

	  int nf = faceused.Size();
	  for (int i = 0; i < 8; i++)
	  {
		FindInnerBoxesRec2(box.childs[i], adfront, faceboxes, faceinds, nf);
	  }
	}



	///
	private void SetInnerBoxesRec(GradingBox box)
	{
	  box.flags.isinner = 1;
	  for (int i = 0; i < 8; i++)
	  {
		if (box.childs[i])
		{
	  ClearFlagsRec(box.childs[i]);
		}
	  }
	}

	///
	private void ClearFlagsRec(GradingBox box)
	{
	  box.flags.cutboundary = 0;
	  box.flags.isinner = 0;
	  for (int i = 0; i < 8; i++)
	  {
		if (box.childs[i])
		{
	  ClearFlagsRec(box.childs[i]);
		}
	  }
	}

	///
	private void ConvexifyRec(GradingBox box)
	{
	  Point < 3> center = box.PMid();

	  double size = 2 * box.h2; // box->x2[0] - box->x1[0];
	  double dx = 0.6 * size;

	  double maxh = box.hopt;

	  for (int i = 0; i < 3; i++)
	  {
	  Point < 3> hp = center;
	  hp(i) += dx;
	  maxh = netgen.GlobalMembers.max2(maxh, GetH(hp));
	  hp(i) = center(i) - dx;
	  maxh = netgen.GlobalMembers.max2(maxh, GetH(hp));
	  }

	  if (maxh < 0.95 * box.hopt)
	  {
		SetH(center, maxh);
	  }

	  for (int i = 0; i < 8; i++)
	  {
		if (box.childs[i])
		{
	  ConvexifyRec(box.childs[i]);
		}
	  }
	}

	private static ostream operator << (ostream ost, LocalH loch)
	{
	  for (int i = 0; i < loch.boxes.Size(); i++)
	  {
		ost << "box[" << i << "] = " << *(loch.boxes[i]);
	  }
	  return ost;
	}
  }

}




namespace netgen
{
  private delegate int testinnerDelegate(Point3d p1);


  private delegate int testinnerDelegate(Point < 2> p1);


}
