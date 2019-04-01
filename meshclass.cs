using Point = netgen.Point;
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
/* File:   meshclass.hpp                                                  */
/* Author: Joachim Schoeberl                                              */
/* Date:   20. Nov. 99                                                    */
/**************************************************************************/

/*
  The mesh class
*/

namespace netgen
{
  public enum resthtype
  {
	  RESTRICTH_FACE,
	  RESTRICTH_EDGE,
		   RESTRICTH_SURFACEELEMENT,
		   RESTRICTH_POINT,
		   RESTRICTH_SEGMENT
  }

//C++ TO C# CONVERTER NOTE: C# has no need of forward class declarations:
//  class HPRefElement;


  /// 2d/3d mesh
  public class Mesh : System.IDisposable
  {
	public typedef.netgen.T_POINTS T_POINTS = new typedef.netgen.T_POINTS();
	// typedef Array<Element2d, 0, SurfaceElementIndex> T_SURFELEMENTS;

	/// point coordinates
	private global::netgen.T_POINTS points = new global::netgen.T_POINTS();

	// The communicator for this mesh. Just a dummy if compiled without MPI.  
	private NgMPI_Comm comm = new NgMPI_Comm();

	/// line-segments at edges
	private Array<Segment, 0, uint> segments = new Array<Segment, 0, uint>();
	/// surface elements, 2d-inner elements
	private Array<Element2d, 0, uint> surfelements = new Array<Element2d, 0, uint>();
	/// volume elements
	private Array<Element, 0, uint> volelements = new Array<Element, 0, uint>();
	/// points will be fixed forever
	private Array<PointIndex> lockedpoints = new Array<PointIndex>();


	/// surface indices at boundary nodes
	// TABLE<int,PointIndex::BASE> surfacesonnode;
	/// boundary edges  (1..normal bedge, 2..segment)
	private INDEX_2_CLOSED_HASHTABLE<int> boundaryedges;
	///
	private INDEX_2_CLOSED_HASHTABLE<int> segmentht;
	///
	private INDEX_3_CLOSED_HASHTABLE<int> surfelementht;

	/// faces of rest-solid
	private Array<Element2d> openelements = new Array<Element2d>();
	/// open segmenets for surface meshing  
	private Array<Segment> opensegments = new Array<Segment>();



	/**
	   Representation of local mesh-size h
	*/
	private LocalH lochfunc;
	///
	private double hglob;
	///
	private double hmin;
	///
	private Array<double> maxhdomain = new Array<double>();

	/**
	   the face-index of the surface element maps into
	   this table.
	*/
	private Array<FaceDescriptor> facedecoding = new Array<FaceDescriptor>();


	/**
	   the edge-index of the line element maps into
	   this table.
	*/
	private Array<EdgeDescriptor> edgedecoding = new Array<EdgeDescriptor>();

	/// sub-domain materials 
	private Array<string> materials = new Array<string>();

	/// labels for boundary conditions
	private Array<string> bcnames = new Array<string>();

	/// labels for co dim 2 bboundary conditions
	private Array<string> cd2names = new Array<string>();

	/// labels for co dim 3 bbboundary conditions
	private Array<string> cd3names = new Array<string>();

	/// Periodic surface, close surface, etc. identifications
	private Identifications ident;


	/// number of vertices (if < 0, use np)
	private int numvertices;

	/// geometric search tree for interval intersection search
	private BoxTree < 3> * elementsearchtree;
	/// time stamp for tree
	private int elementsearchtreets;

	/// element -> face, element -> edge etc ...
	private MeshTopology topology = new MeshTopology();
	/// methods for high order elements
	private CurvedElements curvedelems;

	/// nodes identified by close points 
	private AnisotropicClusters clusters;

	/// space dimension (2 or 3)
	private int dimension;

	/// changed by every minor modification (addpoint, ...)
	private int timestamp;
	/// changed after finishing global algorithm (improve, ...)
	private int majortimestamp;

	/// mesh access semaphors.
	private NgMutex mutex = new NgMutex();
	/// mesh access semaphors.
	private NgMutex majormutex = new NgMutex();

	private SymbolTable< Array<int> > userdata_int = new SymbolTable< Array<int> >();
	private SymbolTable< Array<double> > userdata_double = new SymbolTable< Array<double> >();


	private Array< Point3d > pointcurves = new Array< Point3d >();
	private Array<int> pointcurves_startpoint = new Array<int>();
	private Array<double> pointcurves_red = new Array<double>();
	private Array<double> pointcurves_green = new Array<double>();
	private Array<double> pointcurves_blue = new Array<double>();


	/// start element for point search (GetElementOfPoint)
	private int ps_startelement;


#if PARALLEL
	/// connection to parallel meshes
	private ParallelMeshTopology paralleltop;

#endif


	private NetgenGeometry geometry;

	private void BuildBoundaryEdges()
	{
	  boundaryedges = null;

	  boundaryedges = new INDEX_2_CLOSED_HASHTABLE<int> (3 * (GetNSE() + GetNOpenElements()) + GetNSeg() + 1);


	  for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
	  {
		  Element2d sel = surfelements[sei];
		  if (sel.IsDeleted())
		  {
			  continue;
		  }

		  // int si = sel.GetIndex();

		  if (sel.GetNP() <= 4)
		  {
			for (int j = 0; j < sel.GetNP(); j++)
			{
				INDEX_2 i2 = new INDEX_2();
				i2.I1() = sel.PNumMod(j + 1);
				i2.I2() = sel.PNumMod(j + 2);
				i2.Sort();
				boundaryedges.Set(i2, 1);
			}
		  }
		  else if (sel.GetType() == ELEMENT_TYPE.TRIG6)
		  {
			  for (int j = 0; j < 3; j++)
			  {
				  INDEX_2 i2 = new INDEX_2();
				  i2.I1() = sel[j];
				  i2.I2() = sel[(j + 1) % 3];
				  i2.Sort();
				  boundaryedges.Set(i2, 1);
			  }
		  }
		  else
		  {
			cerr << "illegal element for buildboundaryedges" << "\n";
		  }
	  }


	  for (int i = 0; i < openelements.Size(); i++)
	  {
		  Element2d sel = openelements[i];
		  for (int j = 0; j < sel.GetNP(); j++)
		  {
			  INDEX_2 i2 = new INDEX_2();
			  i2.I1() = sel.PNumMod(j + 1);
			  i2.I2() = sel.PNumMod(j + 2);
			  i2.Sort();
			  boundaryedges.Set(i2, 1);

			  points[sel[j]].SetType(POINTTYPE.FIXEDPOINT);
		  }
	  }

	  for (int i = 0; i < GetNSeg(); i++)
	  {
		  Segment seg = segments[i];
		  INDEX_2 i2 = new INDEX_2(seg[0], seg[1]);
		  i2.Sort();

		  boundaryedges.Set(i2, 2);
		  //segmentht -> Set (i2, i);
	  }


	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool PointContainedIn2DElement(const Point3d & p, double lami[3], const int element, bool consider3D = false) const
	public bool PointContainedIn2DElement(Point3d p, double[] lami, int element, bool consider3D = false)
	{
	  Vec3d col1 = new Vec3d();
	  Vec3d col2 = new Vec3d();
	  Vec3d col3 = new Vec3d();
	  Vec3d rhs = new Vec3d();
	  Vec3d sol = new Vec3d();
	  const double eps = 1e-6;

	  Array<Element2d> loctrigs = new Array<Element2d>();


	  //SZ
	  if (SurfaceElement(element).GetType() == ELEMENT_TYPE.QUAD)
	  {
		  Element2d el = SurfaceElement(element);

		  Point3d p1 = new Point(el.PNum(1));
		  Point3d p2 = new Point(el.PNum(2));
		  Point3d p3 = new Point(el.PNum(3));
		  Point3d p4 = new Point(el.PNum(4));

		  // Coefficients of Bilinear Mapping from Ref-Elem to global Elem
		  // X = a + b x + c y + d x y
		  Vec3d a = p1;
		  Vec3d b = p2 - a;
		  Vec3d c = p4 - a;
		  Vec3d d = p3 - a - b - c;

		  /*cout << "p = " << p << endl;
		  cout << "p1 = " << p1 << endl;
		  cout << "p2 = " << p2 << endl;
		  cout << "p3 = " << p3 << endl;
		  cout << "p4 = " << p4 << endl;
  
		  cout << "a = " << a << endl;
		  cout << "b = " << b << endl;
		  cout << "c = " << c << endl;
		  cout << "d = " << d << endl;*/


		  Vec3d pa = p - a;
		  double dxb = d.X() * b.Y() - d.Y() * b.X();
		  double dxc = d.X() * c.Y() - d.Y() * c.X();
		  double bxc = b.X() * c.Y() - b.Y() * c.X();
		  double bxpa = b.X() * pa.Y() - b.Y() * pa.X();
		  double cxpa = c.X() * pa.Y() - c.Y() * pa.X();
		  double dxpa = d.X() * pa.Y() - d.Y() * pa.X();

		  /*cout << "dxb = " << dxb << endl;
		  cout << "dxc = " << dxc << endl;
		  cout << "bxc = " << bxc << endl;
		  cout << "bxpa = " << bxpa << endl;
		  cout << "cxpa = " << cxpa << endl;
		  cout << "dxpa = " << dxpa << endl;*/

		  /*
		    P = a + b x + c y + d x y
		    1) P1 = a1 + b1 x + c1 y + d1 x y
		    2) P2 = a2 + b2 x + c2 y + d2 x y
		  
		    -> det(x,d) = det(a,d) + det(b,d) x + det(c,d) y
		      -> x = 1/det(b,d) *( det(P-a,d)-det(c,d) y )
		      -> y = 1/det(c,d) *( det(P-a,d)-det(b,d) x )
		  
		    -> x = (P1 - a1 - c1 y)/(b1 + d1 y)
		      -> det(c,d) y**2 + [det(d,P-a) + det(c,b)] y + det(b,P-a) = 0
		    ( same if we express x = (P2 - a2 - c2 y)/(b2 + d2 y) )
  
		    -> y = (P1 - a1 - b1 x)/(c1 + d1 x)
		      -> det(b,d) x**2 + [det(d,P-a) + det(b,c)] x + det(c,P-a) = 0
		    ( same if we express y = (P2 - a2 - b2 x)/(c2 + d2 x)
		   */

		  lami[2] = 0.0;
		  double eps = 1.E-12;
		  double c1;
		  double c2;
		  double r;

		  //First check if point is "exactly" a vertex point
		  Vec3d d1 = p - p1;
		  Vec3d d2 = p - p2;
		  Vec3d d3 = p - p3;
		  Vec3d d4 = p - p4;

		  //cout << " d1 = " << d1 << ", d2 = " << d2 << ", d3 = " << d3 << ", d4 = " << d4 << endl;

		  if (d1.Length2() < netgen.GlobalMembers.sqr(eps) * d2.Length2() && d1.Length2() < netgen.GlobalMembers.sqr(eps) * d3.Length2() && d1.Length2() < netgen.GlobalMembers.sqr(eps) * d4.Length2())
		  {
			  lami[0] = lami[1] = 0.0;
			  return true;
		  }
		  else if (d2.Length2() < netgen.GlobalMembers.sqr(eps) * d1.Length2() && d2.Length2() < netgen.GlobalMembers.sqr(eps) * d3.Length2() && d2.Length2() < netgen.GlobalMembers.sqr(eps) * d4.Length2())
		  {
			  lami[0] = 1.0;
			  lami[1] = 0.0;
			  return true;
		  }
		  else if (d3.Length2() < netgen.GlobalMembers.sqr(eps) * d1.Length2() && d3.Length2() < netgen.GlobalMembers.sqr(eps) * d2.Length2() && d3.Length2() < netgen.GlobalMembers.sqr(eps) * d4.Length2())
		  {
			  lami[0] = lami[1] = 1.0;
			  return true;
		  }
		  else if (d4.Length2() < netgen.GlobalMembers.sqr(eps) * d1.Length2() && d4.Length2() < netgen.GlobalMembers.sqr(eps) * d2.Length2() && d4.Length2() < netgen.GlobalMembers.sqr(eps) * d3.Length2())
		  {
			  lami[0] = 0.0;
			  lami[1] = 1.0;
			  return true;
		  } //if d is nearly 0: solve resulting linear system
		  else if (d.Length2() < netgen.GlobalMembers.sqr(eps) * b.Length2() && d.Length2() < netgen.GlobalMembers.sqr(eps) * c.Length2())
		  {
			  Vec2d sol = new Vec2d();
			  netgen.GlobalMembers.SolveLinearSystemLS(b, c, p - a, sol);
			  lami[0] = sol.X();
			  lami[1] = sol.Y();
		  return netgen.GlobalMembers.ValidBarCoord(lami, eps);
		  } // if dxc is nearly 0: solve resulting linear equation for y and compute x
		  else if (ngsimd.GlobalMembers.fabs(dxc) < netgen.GlobalMembers.sqr(eps))
		  {
			  lami[1] = -bxpa / (dxpa - bxc);
			  lami[0] = (dxpa - dxc * lami[1]) / dxb;
			  return netgen.GlobalMembers.ValidBarCoord(lami, eps);
		  } // if dxb is nearly 0: solve resulting linear equation for x and compute y
		  else if (ngsimd.GlobalMembers.fabs(dxb) < netgen.GlobalMembers.sqr(eps))
		  {
			  lami[0] = -cxpa / (dxpa + bxc);
			  lami[1] = (dxpa - dxb * lami[0]) / dxc;
			  return netgen.GlobalMembers.ValidBarCoord(lami, eps);
		  } //if dxb >= dxc: solve quadratic equation in y and compute x
		  else if (ngsimd.GlobalMembers.fabs(dxb) >= ngsimd.GlobalMembers.fabs(dxc))
		  {
			  c1 = (bxc - dxpa) / dxc;
			  c2 = -bxpa / dxc;
			  r = c1 * c1 / 4.0 - c2;

			  //quadratic equation has only 1 (unstable) solution
			  if (ngsimd.GlobalMembers.fabs(r) < eps) //not eps^2!
			  {
				  lami[1] = -c1 / 2;
				  lami[0] = (dxpa - dxc * lami[1]) / dxb;
				  return netgen.GlobalMembers.ValidBarCoord(lami, eps);
			  }
			  if (r < 0)
			  {
				  return false;
			  }

			  lami[1] = -c1 / 2 + ngsimd.GlobalMembers.sqrt(r);
			  lami[0] = (dxpa - dxc * lami[1]) / dxb;

			  if (netgen.GlobalMembers.ValidBarCoord(lami, eps))
			  {
				  return true;
			  }
			  else
			  {
				  lami[1] = -c1 / 2 - ngsimd.GlobalMembers.sqrt(r);
				  lami[0] = (dxpa - dxc * lami[1]) / dxb;
				  return netgen.GlobalMembers.ValidBarCoord(lami, eps);
			  }
		  } //if dxc > dxb: solve quadratic equation in x and compute y
		  else
		  {
			  c1 = (-bxc - dxpa) / dxb;
			  c2 = -cxpa / dxb;
			  r = c1 * c1 / 4.0 - c2;

			  //quadratic equation has only 1 (unstable) solution
			  if (ngsimd.GlobalMembers.fabs(r) < eps) //not eps^2!
			  {
				  lami[0] = -c1 / 2;
				  lami[1] = (dxpa - dxb * lami[0]) / dxc;
				  return netgen.GlobalMembers.ValidBarCoord(lami, eps);
			  }
			  if (r < 0)
			  {
				  return false;
			  }

			  lami[0] = -c1 / 2 + ngsimd.GlobalMembers.sqrt(r);
			  lami[1] = (dxpa - dxb * lami[0]) / dxc;

			  if (netgen.GlobalMembers.ValidBarCoord(lami, eps))
			  {
				  return true;
			  }
			  else
			  {
				  lami[0] = -c1 / 2 - ngsimd.GlobalMembers.sqrt(r);
				  lami[1] = (dxpa - dxb * lami[0]) / dxc;
				  return netgen.GlobalMembers.ValidBarCoord(lami, eps);
			  }
		  }

		  /*
		  double dxa = d.X()*a.Y()-d.Y()*a.X();
		  double dxp = d.X()*p.Y()-d.Y()*p.X();
  
  
		  double c0,c1,c2; // ,rt;
		  
  
	  Vec3d dp13 = p3-p1;
	  Vec3d dp24 = p4-p2;
	  double d1 = dp13.Length2();
	  double d2 = dp24.Length2();
  
	  // if(fabs(d.X()) <= eps && fabs(d.Y())<= eps)
	  //if (d.Length2() < sqr(eps))
		  if (d.Length2() < sqr(eps)*d1 && d.Length2() < sqr(eps)*d2)
		    {
		  //Solve Linear System
		  Vec2d sol;
		      SolveLinearSystemLS (b, c, p-a, sol);
		      lami[0] = sol.X();
		      lami[1] = sol.Y();
  
		  if(lami[1]<=1.+eps && lami[1]>=0.-eps && lami[0]<=1.+eps && lami[0]>=0.-eps)
			return true;
  
		  
		        //lami[0]=(c.Y()*(p.X()-a.X())-c.X()*(p.Y()-a.Y()))/
		        //(b.X()*c.Y() -b.Y()*c.X());
		      //lami[1]=(-b.Y()*(p.X()-a.X())+b.X()*(p.Y()-a.Y()))/
		       // (b.X()*c.Y() -b.Y()*c.X());
		  
		    }
		  else
		    if(fabs(dxb) <= eps*fabs(dxc))
		      {
			lami[1] = (dxp-dxa)/dxc;
		        if(fabs(b.X()+d.X()*lami[1])>=fabs(b.Y()+d.Y()*lami[1]))
		          lami[0] = (p.X()-a.X() - c.X()*lami[1])/(b.X()+d.X()*lami[1]);
		        else
		          lami[0] = (p.Y()-a.Y() - c.Y()*lami[1])/(b.Y()+d.Y()*lami[1]);
  
			if(lami[1]<=1.+eps && lami[1]>=0.-eps && lami[0]<=1.+eps && lami[0]>=0.-eps)
		  return true;
		      }
		    else
		      if(fabs(dxc) <= eps*fabs(dxb))
		        {
		  lami[0] = (dxp-dxa)/dxb;
		          if(fabs(c.X()+d.X()*lami[0])>=fabs(c.Y()+d.Y()*lami[0]))
		            lami[1] = (p.X()-a.X() - b.X()*lami[0])/(c.X()+d.X()*lami[0]);
		          else
		            lami[1] = (p.Y()-a.Y() - b.Y()*lami[0])/(c.Y()+d.Y()*lami[0]);
  
		  if(lami[1]<=1.+eps && lami[1]>=0.-eps && lami[0]<=1.+eps && lami[0]>=0.-eps)
			return true;
		        }
		      else //Solve quadratic equation
		        {
		  c2 = -d.X()*dxb;
		  c1 = b.X()*dxc - c.X()*dxb + d.X()*(dxp-dxa);
		  c0 = c.X()*(dxp-dxa) + (a.X()-p.X())*dxc;
		  double rt =  c1*c1 - 4*c2*c0;
  
		  if (rt < 0.) return false;
		  lami[1] = (-c1 + sqrt(rt))/2/c2;
  
  
		  if(lami[1]<=1.+eps && lami[1]>=0.-eps)
			{
			  lami[0] = (dxp - dxa -dxb*lami[1])/dxc;
  
			  if(lami[0]<=1.+eps && lami[0]>=0.-eps)
				return true;
			}
		  lami[1] = (-c1 - sqrt(rt))/2/c2;
  
		  lami[0] = (dxp - dxa -dxb*lami[1])/dxc;
  
		  if(lami[1]<=1.+eps && lami[1]>=0.-eps && lami[0]<=1.+eps && lami[0]>=0.-eps)
			return true;
  
		  c2 = d.Y()*dxb;
		  c1 = b.Y()*dxc - c.Y()*dxb + d.Y()*(dxp-dxa);
		  c0 = c.Y()*(dxp -dxa) + (a.Y()-p.Y())*dxc;
		  rt =  c1*c1 - 4*c2*c0;
  
		  if (rt < 0.) return false;
		  lami[1] = (-c1 + sqrt(rt))/2/c2;
  
		  if(lami[1]<=1.+eps && lami[1]>=0.-eps)
			{
			  lami[0] = (dxp - dxa -dxb*lami[1])/dxc;
  
			  if(lami[0]<=1.+eps && lami[0]>=0.-eps)
				return true;
			}
		  lami[1] = (-c1 - sqrt(rt))/2/c2;
  
		  lami[0] = (dxp - dxa -dxb*lami[1])/dxc;
  
		  if(lami[1]<=1.+eps && lami[1]>=0.-eps && lami[0]<=1.+eps && lami[0]>=0.-eps)
			return true;
  
		  c2 = -d.X()*dxc;
		  c1 = -b.X()*dxc + c.X()*dxb + d.X()*(dxp-dxa);
		  c0 = b.X()*(dxp -dxa) + (a.X()-p.X())*dxb;
		  rt =  c1*c1 - 4*c2*c0;
  
		  if (rt < 0.) return false;
		  lami[1] = (-c1 + sqrt(rt))/2/c2;
  
		  if(lami[1]<=1.+eps && lami[1]>=0.-eps)
			{
			  lami[0] = (dxp - dxa -dxc*lami[1])/dxb;
  
			  if(lami[0]<=1.+eps && lami[0]>=0.-eps)
				return true;
			}
		  lami[1] = (-c1 - sqrt(rt))/2/c2;
  
		  lami[0] = (dxp - dxa -dxc*lami[1])/dxb;
  
		  if(lami[1]<=1.+eps && lami[1]>=0.-eps && lami[0]<=1.+eps && lami[0]>=0.-eps)
			return true;
		            }*/


		  //cout << "lam0,1 = " << lami[0] << ", " << lami[1] << endl;

		  /*if( lami[0] <= 1.+eps  && lami[0] >= -eps && lami[1]<=1.+eps && lami[1]>=-eps)
		    {
		      if(consider3D)
		        {
		          Vec3d n = Cross(b,c);
		          lami[2] = 0;
		          for(int i=1; i<=3; i++)
		            lami[2] +=(p.X(i)-a.X(i)-lami[0]*b.X(i)-lami[1]*c.X(i)) * n.X(i);
		          if(lami[2] >= -eps && lami[2] <= eps)
		            return true;
		        }
		      else
		        return true;
			}*/

		  return false;

	  }
	  else
	  {
		  //	  SurfaceElement(element).GetTets (loctets);
		  loctrigs.SetSize(1);
		  loctrigs.Elem(1) = SurfaceElement(element);



		  for (int j = 1; j <= loctrigs.Size(); j++)
		  {
			  Element2d el = loctrigs.Get(j);


			  Point3d p1 = new Point(el.PNum(1));
			  Point3d p2 = new Point(el.PNum(2));
			  Point3d p3 = new Point(el.PNum(3));
			  /*
			    Box3d box;
			    box.SetPoint (p1);
			    box.AddPoint (p2);
			    box.AddPoint (p3);
			    box.AddPoint (p4);
			    if (!box.IsIn (p))
			    continue;
			  */
			  col1 = p2 - p1;
			  col2 = p3 - p1;
			  col3 = netgen.GlobalMembers.Cross(col1,col2);
			  //col3 = Vec3d(0, 0, 1);
			  rhs = p - p1;

			  // int retval =
			  SolveLinearSystem(col1, col2, col3, rhs, sol);

			  //(*testout) << "retval " << retval << endl;

			  //(*testout) << "col1 " << col1 << " col2 " << col2 << " col3 " << col3 << " rhs " << rhs << endl;
			  //(*testout) << "sol " << sol << endl;

			  if (sol.X() >= -eps && sol.Y() >= -eps && sol.X() + sol.Y() <= 1 + eps)
			  {
				  if (!consider3D || (sol.Z() >= -eps && sol.Z() <= eps))
				  {
					  lami[0] = sol.X();
					  lami[1] = sol.Y();
					  lami[2] = sol.Z();

					  return true;
				  }
			  }
		  }
	  }

	  return false;

	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool PointContainedIn3DElement(const Point3d & p, double lami[3], const int element) const
	public bool PointContainedIn3DElement(Point3d p, double[] lami, int element)
	{
	  //bool oldresult = PointContainedIn3DElementOld(p,lami,element);
	  //(*testout) << "old result: " << oldresult
	  //       << " lam " << lami[0] << " " << lami[1] << " " << lami[2] << endl;

	  //if(!curvedelems->IsElementCurved(element-1))
	  //  return PointContainedIn3DElementOld(p,lami,element);


	  const double eps = 1.e-4;
	  Element el = VolumeElement(element);

	  netgen.Point < 3> lam = 0.0;

	  if (el.GetType() == ELEMENT_TYPE.TET || el.GetType() == ELEMENT_TYPE.TET10)
	  {
		  lam = 0.25;
	  }
	  else if (el.GetType() == ELEMENT_TYPE.PRISM)
	  {
		  lam(0) = 0.33;
		  lam(1) = 0.33;
		  lam(2) = 0.5;
	  }
	  else if (el.GetType() == ELEMENT_TYPE.PYRAMID)
	  {
		  lam(0) = 0.4;
		  lam(1) = 0.4;
		  lam(2) = 0.2;
	  }
	  else if (el.GetType() == ELEMENT_TYPE.HEX)
	  {
		  lam = 0.5;
	  }


	  Vec < 3> deltalam,rhs;
	  netgen.Point < 3> x;
	  Mat < 3,3> Jac,Jact;

	  double delta = 1;

	  bool retval;

	  int i = 0;

	  const int maxits = 30;
	  while (delta > 1e-16 && i < maxits)
	  {
		  curvedelems.CalcElementTransformation(lam,element - 1,x,Jac);
		  rhs = p - x;
		  Jac.Solve(rhs,deltalam);

		  lam += deltalam;

		  delta = deltalam.Length2();

		  i++;
		  //(*testout) << "pcie i " << i << " delta " << delta << " p " << p << " x " << x << " lam " << lam << endl;
		  //<< "Jac " << Jac << endl;
	  }

	  if (i == maxits)
	  {
		return false;
	  }

	  for (i = 0; i < 3; i++)
	  {
		lami[i] = lam(i);
	  }



	  if (el.GetType() == ELEMENT_TYPE.TET || el.GetType() == ELEMENT_TYPE.TET10)
	  {
		  retval = (lam(0) > -eps && lam(1) > -eps && lam(2) > -eps && lam(0) + lam(1) + lam(2) < 1 + eps);
	  }
	  else if (el.GetType() == ELEMENT_TYPE.PRISM || el.GetType() == ELEMENT_TYPE.PRISM15)
	  {
		  retval = (lam(0) > -eps && lam(1) > -eps && lam(2) > -eps && lam(2) < 1 + eps && lam(0) + lam(1) < 1 + eps);
	  }
	  else if (el.GetType() == ELEMENT_TYPE.PYRAMID || el.GetType() == ELEMENT_TYPE.PYRAMID13)
	  {
		  retval = (lam(0) > -eps && lam(1) > -eps && lam(2) > -eps && lam(0) + lam(2) < 1 + eps && lam(1) + lam(2) < 1 + eps);
	  }
	  else if (el.GetType() == ELEMENT_TYPE.HEX || el.GetType() == ELEMENT_TYPE.HEX20)
	  {
		  retval = (lam(0) > -eps && lam(0) < 1 + eps && lam(1) > -eps && lam(1) < 1 + eps && lam(2) > -eps && lam(2) < 1 + eps);
	  }
	  else
	  {
		throw new Exception("Da haun i wos vagessn");
	  }

	  return retval;
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool PointContainedIn3DElementOld(const Point3d & p, double lami[3], const int element) const
	public bool PointContainedIn3DElementOld(Point3d p, double[] lami, int element)
	{
	  Vec3d col1 = new Vec3d();
	  Vec3d col2 = new Vec3d();
	  Vec3d col3 = new Vec3d();
	  Vec3d rhs = new Vec3d();
	  Vec3d sol = new Vec3d();
	  const double eps = 1.e-4;

	  Array<Element> loctets = new Array<Element>();

	  VolumeElement(element).GetTets(loctets);

	  for (int j = 1; j <= loctets.Size(); j++)
	  {
		  Element el = loctets.Get(j);

		  Point3d p1 = new Point(el.PNum(1));
		  Point3d p2 = new Point(el.PNum(2));
		  Point3d p3 = new Point(el.PNum(3));
		  Point3d p4 = new Point(el.PNum(4));

		  Box3d box = new Box3d();
		  box.SetPoint(p1);
		  box.AddPoint(p2);
		  box.AddPoint(p3);
		  box.AddPoint(p4);
		  if (box.IsIn(p) == 0)
		  {
			continue;
		  }

		  col1 = p2 - p1;
		  col2 = p3 - p1;
		  col3 = p4 - p1;
		  rhs = p - p1;

		  SolveLinearSystem(col1, col2, col3, rhs, sol);

		  if (sol.X() >= -eps && sol.Y() >= -eps && sol.Z() >= -eps && sol.X() + sol.Y() + sol.Z() <= 1 + eps)
		  {
			  Array<Element> loctetsloc = new Array<Element>();
			  Array<netgen.Point < 3>> pointsloc = new Array<netgen.Point < 3>>();

			  VolumeElement(element).GetTetsLocal(loctetsloc);
			  VolumeElement(element).GetNodesLocalNew(pointsloc);

			  Element le = loctetsloc.Get(j);


			  Point3d pp = pointsloc.Get(le.PNum(1)) + sol.X() * new Vec3d(pointsloc.Get(le.PNum(1)), pointsloc.Get(le.PNum(2))) + sol.Y() * new Vec3d(pointsloc.Get(le.PNum(1)), pointsloc.Get(le.PNum(3))) + sol.Z() * new Vec3d(pointsloc.Get(le.PNum(1)), pointsloc.Get(le.PNum(4)));

			  lami[0] = pp.X();
			  lami[1] = pp.Y();
			  lami[2] = pp.Z();
			  return true;
		  }
	  }
	  return false;
	}


	// store coarse mesh before hp-refinement
	public Array<HPRefElement>[] hpelements;
	public Mesh coarsemesh;


	/// number of refinement levels
	public int mglevels;
	/// refinement hierarchy
	public Array<PointIndices < 2>,PointIndex.BASE> mlbetweennodes = new Array<PointIndices < 2>,PointIndex.BASE>();
	/// parent element of volume element
	public Array<int> mlparentelement = new Array<int>();
	/// parent element of surface element
	public Array<int> mlparentsurfaceelement = new Array<int>();



	///
	public Mesh()
	{
		this.topology = new netgen.MeshTopology(this);
		this.surfarea = new netgen.Mesh.CSurfaceArea(this);
	  // volelements.SetName ("vol elements");
	  // surfelements.SetName ("surf elements");
	  // points.SetName ("meshpoints");

	  boundaryedges = null;
	  surfelementht = null;
	  segmentht = null;

	  lochfunc = null;
	  mglevels = 1;
	  elementsearchtree = null;
	  elementsearchtreets = netgen.GlobalMembers.NextTimeStamp();
	  majortimestamp = timestamp = netgen.GlobalMembers.NextTimeStamp();
	  hglob = 1e10;
	  hmin = 0;
	  numvertices = -1;
	  dimension = 3;

	  // topology = new MeshTopology (*this);
	  curvedelems = new CurvedElements(this);
	  clusters = new AnisotropicClusters(this);
	  ident = new Identifications(this);

	  hpelements = null;
	  coarsemesh = null;

	  ps_startelement = 0;

	  geomtype = GEOM_TYPE.NO_GEOM;

	  bcnames.SetSize(0);
	  cd2names.SetSize(0);

	  // this->comm = netgen :: ng_comm;
#if PARALLEL
	  paralleltop = new ParallelMeshTopology(this);
#endif
	}

	///
	public void Dispose()
	{
	  // cout << "******************** deleting Mesh **********" << endl;
	  if (lochfunc != null)
	  {
		  lochfunc.Dispose();
	  }
	  if (boundaryedges != null)
	  {
		  boundaryedges.Dispose();
	  }
	  if (surfelementht != null)
	  {
		  surfelementht.Dispose();
	  }
	  if (segmentht != null)
	  {
		  segmentht.Dispose();
	  }
	  if (curvedelems != null)
	  {
		  curvedelems.Dispose();
	  }
	  if (clusters != null)
	  {
		  clusters.Dispose();
	  }
	  // delete topology;
	  if (ident != null)
	  {
		  ident.Dispose();
	  }
	  if (elementsearchtree != null)
	  {
		  elementsearchtree.Dispose();
	  }
	  if (coarsemesh != null)
	  {
		  coarsemesh.Dispose();
	  }
	  if (hpelements != null)
	  {
		  hpelements.Dispose();
	  }

	  for (int i = 0; i < materials.Size(); i++)
	  {
		if (materials[i] != null)
		{
			materials[i].Dispose();
		}
	  }
	  for (int i = 0; i < userdata_int.Size(); i++)
	  {
		if (userdata_int[i] != null)
		{
			userdata_int[i].Dispose();
		}
	  }
	  for (int i = 0; i < userdata_double.Size(); i++)
	  {
		if (userdata_double[i] != null)
		{
			userdata_double[i].Dispose();
		}
	  }

	  for (int i = 0; i < bcnames.Size(); i++)
	  {
		if (bcnames[i] != null)
		{
			bcnames[i].Dispose();
		}
	  }

	  for (int i = 0; i < cd2names.Size(); i++)
	  {
		if (cd2names[i] != null)
		{
			cd2names[i].Dispose();
		}
	  }

#if PARALLEL
	  if (paralleltop != null)
	  {
		  paralleltop.Dispose();
	  }
#endif
	}

//C++ TO C# CONVERTER NOTE: This 'CopyFrom' method was converted from the original copy assignment operator:
//ORIGINAL LINE: Mesh & operator = (const Mesh & mesh2)
	public Mesh CopyFrom (Mesh mesh2)
	{
//C++ TO C# CONVERTER TODO TASK: The following line was determined to be a copy assignment (rather than a reference assignment) - this should be verified and a 'CopyFrom' method should be created:
//ORIGINAL LINE: points = mesh2.points;
	  points.CopyFrom(mesh2.points);
	  // eltyps = mesh2.eltyps;
	  segments.CopyFrom(mesh2.segments);
	  surfelements.CopyFrom(mesh2.surfelements);
	  volelements.CopyFrom(mesh2.volelements);
	  lockedpoints.CopyFrom(mesh2.lockedpoints);
	  facedecoding.CopyFrom(mesh2.facedecoding);
	  dimension = mesh2.dimension;

	  bcnames.SetSize(mesh2.bcnames.Size());
	  for (int i = 0; i < mesh2.bcnames.Size(); i++)
	  {
		if (mesh2.bcnames[i])
		{
			bcnames[i] = new string(*mesh2.bcnames[i]);
		}
		else
		{
			bcnames[i] = 0;
		}
	  }

	  cd2names.SetSize(mesh2.cd2names.Size());
	  for (int i = 0; i < mesh2.cd2names.Size(); i++)
	  {
		if (mesh2.cd2names[i])
		{
			cd2names[i] = new string(*mesh2.cd2names[i]);
		}
		else
		{
			cd2names[i] = 0;
		}
	  }

	  return this;
	}

	///
	public void DeleteMesh()
	{
	  NgLock @lock = new NgLock(object);
	  @lock.Lock();
	  points.SetSize(0);
	  segments.SetSize(0);
	  surfelements.SetSize(0);
	  volelements.SetSize(0);
	  lockedpoints.SetSize(0);
	  // surfacesonnode.SetSize(0);

	  boundaryedges = null;
	  boundaryedges = null;

	  openelements.SetSize(0);
	  facedecoding.SetSize(0);

	  ident = null;
	  ident = new Identifications(this);
	  // delete topology;
	  // topology = new MeshTopology (*this);
	  topology = new MeshTopology(this);
	  curvedelems = null;
	  curvedelems = new CurvedElements(this);
	  if (clusters != null)
	  {
		  clusters.Dispose();
	  }
	  clusters = new AnisotropicClusters(this);

	  for (int i = 0; i < bcnames.Size(); i++)
	  {
		if (bcnames[i])
		{
			bcnames[i] = null;
		}
	  }
	  for (int i = 0; i < cd2names.Size(); i++)
	  {
		if (cd2names[i])
		{
			cd2names[i] = null;
		}
	  }

#if PARALLEL
	  if (paralleltop != null)
	  {
		  paralleltop.Dispose();
	  }
	  paralleltop = new ParallelMeshTopology(this);
#endif

	  @lock.UnLock();

	  timestamp = netgen.GlobalMembers.NextTimeStamp();
	}

	///
	public void ClearSurfaceElements()
	{
	  surfelements.SetSize(0);
	  for (int i = 0; i < facedecoding.Size(); i++)
	  {
		facedecoding[i].firstelement = -1;
	  }

	  timestamp = netgen.GlobalMembers.NextTimeStamp();
	}

	///
	public DLL_HEADER void ClearVolumeElements()
	{
	  volelements.SetSize(0);
	  timestamp = netgen.GlobalMembers.NextTimeStamp();
	}

	///
	public DLL_HEADER void ClearSegments()
	{
	  segments.SetSize(0);
	  timestamp = netgen.GlobalMembers.NextTimeStamp();
	}

	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool TestOk() const
	public bool TestOk()
	{
	  for (ElementIndex ei = 0; ei < volelements.Size(); ei++)
	  {
		  for (int j = 0; j < 4; j++)
		  {
			if (this[ei][j] <= PointIndex.BASE-1)
			{
				(*testout) << "El " << ei << " has 0 nodes: ";
				for (int k = 0; k < 4; k++)
				{
				  (*testout) << this[ei][k];
				}
				break;
			}
		  }
	  }
	  CheckMesh3D(this);
	  return true;
	}

	public void SetAllocSize(int nnodes, int nsegs, int nsel, int nel)
	{
	  points.SetAllocSize(nnodes);
	  segments.SetAllocSize(nsegs);
	  surfelements.SetAllocSize(nsel);
	  volelements.SetAllocSize(nel);
	}


	public PointIndex AddPoint(Point3d p, int layer = 1)
	{
	  return new netgen.PointIndex(AddPoint(p, layer, POINTTYPE.INNERPOINT));
	  /*
	  NgLock lock(mutex);
	  lock.Lock();
  
	  timestamp = NextTimeStamp();
  
	  PointIndex pi = points.End();
	  points.Append ( MeshPoint (p, layer, INNERPOINT) );
  
	  lock.UnLock();
  
	  return pi;
	  */
	}

	public PointIndex AddPoint(Point3d p, int layer, POINTTYPE type)
	{
	  NgLock @lock = new NgLock(object);
	  @lock.Lock();

	  timestamp = netgen.GlobalMembers.NextTimeStamp();

	  PointIndex pi = points.End();
	  points.Append(new MeshPoint(p, layer, type));

	  @lock.UnLock();

	  return new netgen.PointIndex(pi);
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: auto GetNP() const
//C++ TO C# CONVERTER TODO TASK: The return type of the following function could not be determined:
	public auto GetNP()
	{
		return points.Size();
	}

	// [[deprecated("Use Point(PointIndex) instead of int !")]]        
	public MeshPoint Point(int i)
	{
		return points.Elem(i);
	}
	public MeshPoint Point(PointIndex pi)
	{
		return points[pi];
	}
	// [[deprecated("Use Point(PointIndex) instead of int !")]]            
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const MeshPoint & Point(int i) const
	public MeshPoint Point(int i)
	{
		return points.Get(i);
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const MeshPoint & Point(PointIndex pi) const
	public MeshPoint Point(PointIndex pi)
	{
		return points[pi];
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const MeshPoint & operator [] (PointIndex pi) const
	public MeshPoint this [PointIndex pi]
	{
		get
		{
			return points[pi];
		}
		set
		{
			points[pi] = value;
		}
	}
	public MeshPoint this [PointIndex pi]
	{
		get
		{
			return points[pi];
		}
		set
		{
			points[pi] = value;
		}
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const::netgen::T_POINTS & Points() const
	public netgen.T_POINTS Points()
	{
		return new global::netgen.T_POINTS(points);
	}
	public global::netgen.T_POINTS Points()
	{
		return new global::netgen.T_POINTS(points);
	}



	/*
  #ifdef PARALLEL
	PointIndex Mesh :: AddPoint (const Point3d & p, bool isghost,  int layer)
	{ 
	  NgLock lock(mutex);
	  lock.Lock();
  
	  timestamp = NextTimeStamp();
  
	  PointIndex pi = points.Size() + PointIndex::BASE;
	  points.Append ( MeshPoint (p, layer, INNERPOINT) ); 
  
	  lock.UnLock();
  
	  return pi;
	}
  
	PointIndex Mesh :: AddPoint (const Point3d & p, bool isghost, int layer, POINTTYPE type)
	{ 
	  NgLock lock(mutex);
	  lock.Lock();
  
	  timestamp = NextTimeStamp();
  
	  PointIndex pi = points.Size() + PointIndex::BASE;
	  points.Append ( MeshPoint (p, layer, type) ); 
  
	  lock.UnLock();
  
	  return pi;
	}
  
  #endif
	*/


	public SegmentIndex AddSegment(Segment s)
	{
	  NgLock @lock = new NgLock(object);
	  @lock.Lock();
	  timestamp = netgen.GlobalMembers.NextTimeStamp();

	  int maxn = netgen.GlobalMembers.max2(new netgen.Segment(s[0]), new netgen.Segment(s[1]));
	  maxn += 1 - PointIndex.BASE;

	  /*
	    if (maxn > ptyps.Size())
	    {
	    int maxo = ptyps.Size();
	    ptyps.SetSize (maxn);
	    for (int i = maxo; i < maxn; i++)
	    ptyps[i] = INNERPOINT;
	    }
  
	    if (ptyps[s[0]] > EDGEPOINT) ptyps[s[0]] = EDGEPOINT;
	    if (ptyps[s[1]] > EDGEPOINT) ptyps[s[1]] = EDGEPOINT;
	  */

	  if (maxn <= points.Size())
	  {
		  if (points[s[0]].Type() > POINTTYPE.EDGEPOINT)
		  {
			points[s[0]].SetType(POINTTYPE.EDGEPOINT);
		  }
		  if (points[s[1]].Type() > POINTTYPE.EDGEPOINT)
		  {
			points[s[1]].SetType(POINTTYPE.EDGEPOINT);
		  }
	  }
	  /*
	    else
	    {
	    cerr << "edge points nrs > points.Size" << endl;
	    }
	  */

	  SegmentIndex si = segments.Size();
	  segments.Append(s);

	  @lock.UnLock();
//C++ TO C# CONVERTER TODO TASK: The following line was determined to contain a copy constructor call - this should be verified and a copy constructor should be created:
//ORIGINAL LINE: return si;
	  return new netgen.SegmentIndex(si);
	}

	public void DeleteSegment(int segnr)
	{
	  segments.Elem(segnr)[0].Invalidate();
	  segments.Elem(segnr)[1].Invalidate();
	}
	/*
	void FullDeleteSegment (int segnr)  // von wem ist das ???
	{
	  segments.Delete(segnr-PointIndex::BASE);
	}
	*/

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetNSeg() const
	public int GetNSeg()
	{
		return segments.Size();
	}
	// [[deprecated("Use LineSegment(SegmentIndex) instead of int !")]]                
	public Segment LineSegment(int i)
	{
		return segments.Elem(i);
	}
	// [[deprecated("Use LineSegment(SegmentIndex) instead of int !")]]                    
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const Segment & LineSegment(int i) const
	public Segment LineSegment(int i)
	{
		return segments.Get(i);
	}

	public Segment LineSegment(SegmentIndex si)
	{
		return segments[si];
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const Segment & LineSegment(SegmentIndex si) const
	public Segment LineSegment(SegmentIndex si)
	{
		return segments[si];
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const Segment & operator [] (SegmentIndex si) const
	public Segment this [SegmentIndex si]
	{
		get
		{
			return segments[si];
		}
		set
		{
			segments[si] = value;
		}
	}
	public Segment this [SegmentIndex si]
	{
		get
		{
			return segments[si];
		}
		set
		{
			segments[si] = value;
		}
	}

	/*
	const Array<Segment> & LineSegments() const { return segments; }
	Array<Segment> & LineSegments() { return segments; }
	*/
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const auto & LineSegments() const
	public Array<Segment, 0, uint> LineSegments()
	{
		return new Array<Segment, 0, uint>(segments);
	}
	public Array<Segment, 0, uint> LineSegments()
	{
		return new Array<Segment, 0, uint>(segments);
	}

	public Array<Element0d> pointelements = new Array<Element0d>(); // only via python interface

	public SurfaceElementIndex AddSurfaceElement(Element2d el)
	{
	  NgLock @lock = new NgLock(object);
	  @lock.Lock();
	  timestamp = netgen.GlobalMembers.NextTimeStamp();

	  int maxn = el[0];
	  for (int i = 1; i < el.GetNP(); i++)
	  {
		if (el[i] > maxn)
		{
			maxn = el[i];
		}
	  }

	  maxn += 1 - PointIndex.BASE;

	  /*
	    if (maxn > ptyps.Size())
	    {
	    int maxo = ptyps.Size();
	    ptyps.SetSize (maxn);
	    for (i = maxo+PointIndex::BASE;
	    i < maxn+PointIndex::BASE; i++)
	    ptyps[i] = INNERPOINT;
  
	    }
	  */
	  if (maxn <= points.Size())
	  {
		  for (int i = 0; i < el.GetNP(); i++)
		  {
			if (points[el[i]].Type() > POINTTYPE.SURFACEPOINT)
			{
			  points[el[i]].SetType(POINTTYPE.SURFACEPOINT);
			}
		  }
	  }
	  /*
	    else
	    {
	    cerr << "surf points nrs > points.Size" << endl;
	    }
	  */

	  SurfaceElementIndex si = surfelements.Size();
	  surfelements.Append(el);

	  if (el.index <= 0 || el.index > facedecoding.Size())
	  {
		cerr << "has no facedecoding: fd.size = " << facedecoding.Size() << ", ind = " << el.index << "\n";
	  }

	  surfelements.Last().next = facedecoding[el.index - 1].firstelement;
	  facedecoding[el.index - 1].firstelement = si;

	  if (SurfaceArea().Valid())
	  {
		SurfaceArea().Add(el);
	  }

	  @lock.UnLock();
//C++ TO C# CONVERTER TODO TASK: The following line was determined to contain a copy constructor call - this should be verified and a copy constructor should be created:
//ORIGINAL LINE: return si;
	  return new netgen.SurfaceElementIndex(si);
	}

	// write to pre-allocated container, thread-safe
	public void SetSurfaceElement(SurfaceElementIndex sei, Element2d el)
	{
	  int maxn = el[0];
	  for (int i = 1; i < el.GetNP(); i++)
	  {
		if (el[i] > maxn)
		{
			maxn = el[i];
		}
	  }

	  maxn += 1 - PointIndex.BASE;

	  if (maxn <= points.Size())
	  {
		  for (int i = 0; i < el.GetNP(); i++)
		  {
			if (points[el[i]].Type() > POINTTYPE.SURFACEPOINT)
			{
			  points[el[i]].SetType(POINTTYPE.SURFACEPOINT);
			}
		  }
	  }

	  surfelements[sei] = el;
	  if (el.index > facedecoding.Size())
	  {
		cerr << "has no facedecoding: fd.size = " << facedecoding.Size() << ", ind = " << el.index << "\n";
	  }

	  // add lock-free to list ... slow, call RebuildSurfaceElementLists later
	  /*
	  surfelements[sei].next = facedecoding[el.index-1].firstelement;
	  auto & head = reinterpret_cast<atomic<SurfaceElementIndex>&> (facedecoding[el.index-1].firstelement);
	  while (!head.compare_exchange_weak (surfelements[sei].next, sei))
	    ;
	  */

	  /*
	  if (SurfaceArea().Valid())
	    SurfaceArea().Add (el);
	  */
	}

	// [[deprecated("Use DeleteSurfaceElement(SurfaceElementIndex) instead of int !")]]
	public void DeleteSurfaceElement(int eli)
	{
	  surfelements.Elem(eli).Delete();
	  surfelements.Elem(eli).PNum(1).Invalidate();
	  surfelements.Elem(eli).PNum(2).Invalidate();
	  surfelements.Elem(eli).PNum(3).Invalidate();
	  timestamp = netgen.GlobalMembers.NextTimeStamp();
	}

	public void DeleteSurfaceElement(SurfaceElementIndex eli)
	{
	  foreach (var p in surfelements[eli].PNums())
	  {
		  p.Invalidate();
	  }
	  surfelements[eli].Delete();
	  timestamp = netgen.GlobalMembers.NextTimeStamp();
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: auto GetNSE() const
//C++ TO C# CONVERTER TODO TASK: The return type of the following function could not be determined:
	public auto GetNSE()
	{
		return surfelements.Size();
	}

	// [[deprecated("Use SurfaceElement(SurfaceElementIndex) instead of int !")]]    
	public Element2d SurfaceElement(int i)
	{
		return surfelements.Elem(i);
	}
	// [[deprecated("Use SurfaceElement(SurfaceElementIndex) instead of int !")]]        
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const Element2d & SurfaceElement(int i) const
	public Element2d SurfaceElement(int i)
	{
		return surfelements.Get(i);
	}
	public Element2d SurfaceElement(SurfaceElementIndex i)
	{
		return surfelements[i];
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const Element2d & SurfaceElement(SurfaceElementIndex i) const
	public Element2d SurfaceElement(SurfaceElementIndex i)
	{
		return surfelements[i];
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const Element2d & operator [] (SurfaceElementIndex ei) const
	public Element2d this [SurfaceElementIndex ei]
	{
		get
		{
			return surfelements[ei];
		}
		set
		{
			surfelements[ei] = value;
		}
	}
	public Element2d this [SurfaceElementIndex ei]
	{
		get
		{
			return surfelements[ei];
		}
		set
		{
			surfelements[ei] = value;
		}
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const Array<Element2d, 0, uint> & SurfaceElements() const
	public Array<Element2d, 0, uint> SurfaceElements()
	{
		return new Array<Element2d, 0, uint>(surfelements);
	}
	public Array<Element2d, 0, uint> SurfaceElements()
	{
		return new Array<Element2d, 0, uint>(surfelements);
	}


	public void RebuildSurfaceElementLists()
	{
	  for (int i = 0; i < facedecoding.Size(); i++)
	  {
		facedecoding[i].firstelement = -1;
	  }
	  for (int i = surfelements.Size() - 1; i >= 0; i--)
	  {
		  int ind = surfelements[i].GetIndex();
		  surfelements[i].next = facedecoding[ind - 1].firstelement;
		  facedecoding[ind - 1].firstelement = i;
	  }
	}

//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	private int GetSurfaceElementsOfFace_timer = NgProfiler.CreateTimer("GetSurfaceElementsOfFace");

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetSurfaceElementsOfFace(int facenr, Array<SurfaceElementIndex> & sei) const
	public void GetSurfaceElementsOfFace(int facenr, Array<SurfaceElementIndex> sei)
	{
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	  static int timer = NgProfiler::CreateTimer("GetSurfaceElementsOfFace");
	  NgProfiler.RegionTimer reg = new NgProfiler.RegionTimer(GetSurfaceElementsOfFace_timer);

	   /*
	   sei.SetSize (0);
	   for (SurfaceElementIndex i = 0; i < GetNSE(); i++)
	   {
	      if ( (*this)[i].GetIndex () == facenr && (*this)[i][0] >= PointIndex::BASE &&
	         !(*this)[i].IsDeleted() )
	      {
	         sei.Append (i);
	      }
	   }
	   */

	   /* Philippose - 01/10/2009
	   Commented out the following lines, and activated the originally
	   commented out lines above because of a bug which causes corruption
	   of the variable "facedecoding" when a mesh is converted to second order
	   */

	   //      int size1 = sei.Size();
	   sei.SetSize(0);

	   SurfaceElementIndex si = facedecoding[facenr - 1].firstelement;
	   while (si != -1)
	   {
		  if (this[si].GetIndex() == facenr && this[si][0] >= PointIndex.BASE && !this[si].IsDeleted())
		  {
			 sei.Append(si);
		  }

		  si = this[si].next;
	   }

	   /*
	   // *testout << "with list = " << endl << sei << endl;
  
	   if (size1 != sei.Size())
	   {
	      cout << "size mismatch" << endl;
	      exit(1);
	   }
	   */
	}

	public ElementIndex AddVolumeElement(Element el)
	{
	  NgLock @lock = new NgLock(object);
	  @lock.Lock();

	  int maxn = el[0];
	  for (int i = 1; i < el.GetNP(); i++)
	  {
		if (el[i] > maxn)
		{
			maxn = el[i];
		}
	  }

	  maxn += 1 - PointIndex.BASE;

	  /*
	    if (maxn > ptyps.Size())
	    {
	    int maxo = ptyps.Size();
	    ptyps.SetSize (maxn);
	    for (i = maxo+PointIndex::BASE;
	    i < maxn+PointIndex::BASE; i++)
	    ptyps[i] = INNERPOINT;
	    }
	  */
	  /*
	    if (maxn > points.Size())
	    {
	    cerr << "add vol element before point" << endl;
	    }
	  */

	  int ve = volelements.Size();

	  volelements.Append(el);
	  volelements.Last().flags.illegal_valid = 0;

	  // while (volelements.Size() > eltyps.Size())
	  // eltyps.Append (FREEELEMENT);

	  timestamp = netgen.GlobalMembers.NextTimeStamp();

	  @lock.UnLock();
	  return ve;
	}

	// write to pre-allocated container, thread-safe
	public void SetVolumeElement(ElementIndex ei, Element el)
	{
	  /*
	  int maxn = el[0];
	  for (int i = 1; i < el.GetNP(); i++)
	    if (el[i] > maxn) maxn = el[i];
  
	  maxn += 1-PointIndex::BASE;
	  */

	  volelements[ei] = el;
	  volelements.Last().flags.illegal_valid = 0;
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: auto GetNE() const
//C++ TO C# CONVERTER TODO TASK: The return type of the following function could not be determined:
	public auto GetNE()
	{
		return volelements.Size();
	}

	// [[deprecated("Use VolumeElement(ElementIndex) instead of int !")]]    
	public Element VolumeElement(int i)
	{
		return volelements.Elem(i);
	}
	// [[deprecated("Use VolumeElement(ElementIndex) instead of int !")]]        
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const Element & VolumeElement(int i) const
	public Element VolumeElement(int i)
	{
		return volelements.Get(i);
	}
	public Element VolumeElement(ElementIndex i)
	{
		return volelements[i];
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const Element & VolumeElement(ElementIndex i) const
	public Element VolumeElement(ElementIndex i)
	{
		return volelements[i];
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const Element & operator [] (ElementIndex ei) const
	public Element this [ElementIndex ei]
	{
		get
		{
			return volelements[ei];
		}
		set
		{
			volelements[ei] = value;
		}
	}
	public Element this [ElementIndex ei]
	{
		get
		{
			return volelements[ei];
		}
		set
		{
			volelements[ei] = value;
		}
	}



//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: ELEMENTTYPE ElementType(ElementIndex i) const
	public ELEMENTTYPE ElementType(ElementIndex i)
	{
		return (volelements[i].flags.@fixed) ? ELEMENTTYPE.FIXEDELEMENT : ELEMENTTYPE.FREEELEMENT;
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const auto & VolumeElements() const
	public Array<Element, 0, uint> VolumeElements()
	{
		return new Array<Element, 0, uint>(volelements);
	}
	public Array<Element, 0, uint> VolumeElements()
	{
		return new Array<Element, 0, uint>(volelements);
	}

	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: double ElementError(int eli, const MeshingParameters & mp) const
	public double ElementError(int eli, MeshingParameters mp)
	{
	  Element el = volelements.Get(eli);
	  return CalcTetBadness(points.Get(el[0]), points.Get(el[1]), points.Get(el[2]), points.Get(el[3]), -1, mp);
	}

	/// 
	public void AddLockedPoint(PointIndex pi)
	{
	  lockedpoints.Append(pi);
	}

	///
	public void ClearLockedPoints()
	{
	  lockedpoints.SetSize(0);
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const auto & LockedPoints() const
	public Array<PointIndex> LockedPoints()
	{
		return new Array<PointIndex>(lockedpoints);
	}

	/// Returns number of domains
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetNDomains() const
	public int GetNDomains()
	{
	  int ndom = 0;

	  for (int k = 0; k < facedecoding.Size(); k++)
	  {
		  if (facedecoding[k].DomainIn() > ndom)
		  {
			ndom = facedecoding[k].DomainIn();
		  }
		  if (facedecoding[k].DomainOut() > ndom)
		  {
			ndom = facedecoding[k].DomainOut();
		  }
	  }

	  return ndom;
	}

	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetDimension() const
	public int GetDimension()
	{
		return dimension;
	}
	public void SetDimension(int dim)
	{
		dimension = dim;
	}

	/// sets internal tables
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
//C++ TO C# CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
	Timer CalcSurfacesOfNode_t("Mesh::CalcSurfacesOfNode");

	public void CalcSurfacesOfNode()
	{
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	  static Timer t("Mesh::CalcSurfacesOfNode");
	  RegionTimer reg = new RegionTimer(CalcSurfacesOfNode_t);
	  // surfacesonnode.SetSize (GetNP());
	  TABLE<int,PointIndex.BASE> surfacesonnode = new TABLE<int,PointIndex.BASE>(GetNP());

	  boundaryedges = null;
	  boundaryedges = null;

	  surfelementht = null;
	  surfelementht = null;
	  segmentht = null;

	  /*
	    surfelementht = new INDEX_3_HASHTABLE<int> (GetNSE()/4 + 1);
	    segmentht = new INDEX_2_HASHTABLE<int> (GetNSeg() + 1);
	  */

	  if (dimension == 3)
	  {
		surfelementht = new INDEX_3_CLOSED_HASHTABLE<int> (3 * GetNSE() + 1);
	  }
	  segmentht = new INDEX_2_CLOSED_HASHTABLE<int> ((uint)(3 * GetNSeg() + 1));

	  if (dimension == 3)
	  {
	  for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
	  {
		  Element2d sel = surfelements[sei];
		  if (sel.IsDeleted())
		  {
			  continue;
		  }

		  int si = sel.GetIndex();

		  for (int j = 0; j < sel.GetNP(); j++)
		  {
			  PointIndex pi = sel[j];
			  if (!surfacesonnode[pi].Contains(si))
			  {
				surfacesonnode.Add(pi, si);
			  }
			  /*
			  bool found = 0;
			  for (int k = 0; k < surfacesonnode[pi].Size(); k++)
			    if (surfacesonnode[pi][k] == si)
			      {
			        found = 1;
			        break;
			      }
  
			  if (!found)
			    surfacesonnode.Add (pi, si);
			  */
		  }
	  }
	  }
	  /*
	    for (sei = 0; sei < GetNSE(); sei++)
	    {
	    const Element2d & sel = surfelements[sei];
	    if (sel.IsDeleted()) continue;
  
	    INDEX_3 i3;
	    i3.I1() = sel.PNum(1);
	    i3.I2() = sel.PNum(2);
	    i3.I3() = sel.PNum(3);
	    i3.Sort();
	    surfelementht -> PrepareSet (i3);
	    }
  
	    surfelementht -> AllocateElements();
	  */

	  if (dimension == 3)
	  {
	  for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
	  {
		  Element2d sel = surfelements[sei];
		  if (sel.IsDeleted())
		  {
			  continue;
		  }

		  INDEX_3 i3 = new INDEX_3();
		  i3.I1() = sel.PNum(1);
		  i3.I2() = sel.PNum(2);
		  i3.I3() = sel.PNum(3);
		  i3.Sort();
		  surfelementht.Set(i3, sei); // war das wichtig ???    sel.GetIndex());
	  }
	  }

	  // int np = GetNP();

	  if (dimension == 3)
	  {
		  for (PointIndex pi = points.Begin(); pi < points.End(); pi++)
		  {
			points[pi].SetType(POINTTYPE.INNERPOINT);
		  }

		  if (GetNFD() == 0)
		  {
			  for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
			  {
				  Element2d sel = surfelements[sei];
				  if (sel.IsDeleted())
				  {
					  continue;
				  }
				  for (int j = 0; j < sel.GetNP(); j++)
				  {
//C++ TO C# CONVERTER TODO TASK: The following line was determined to contain a copy constructor call - this should be verified and a copy constructor should be created:
//ORIGINAL LINE: PointIndex pi = SurfaceElement(sei)[j];
					  PointIndex pi = SurfaceElement(new netgen.SurfaceElementIndex(sei))[j];
					  points[pi].SetType(POINTTYPE.FIXEDPOINT);
				  }
			  }
		  }
		  else
		  {
			  for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
			  {
				  Element2d sel = surfelements[sei];
				  if (sel.IsDeleted())
				  {
					  continue;
				  }
				  for (int j = 0; j < sel.GetNP(); j++)
				  {
					  PointIndex pi = sel[j];
					  int ns = surfacesonnode[pi].Size();
					  if (ns == 1)
					  {
						points[pi].SetType(POINTTYPE.SURFACEPOINT);
					  }
					  if (ns == 2)
					  {
						points[pi].SetType(POINTTYPE.EDGEPOINT);
					  }
					  if (ns >= 3)
					  {
						points[pi].SetType(POINTTYPE.FIXEDPOINT);
					  }
				  }
			  }
		  }
	  }

	  for (int i = 0; i < segments.Size(); i++)
	  {
		  Segment seg = segments[i];
		  for (int j = 1; j <= 2; j++)
		  {
			  PointIndex hi = (j == 1) ? seg[0] : seg[1];
			  if (points[hi].Type() == POINTTYPE.INNERPOINT || points[hi].Type() == POINTTYPE.SURFACEPOINT)
			  {
				points[hi].SetType(POINTTYPE.EDGEPOINT);
			  }
		  }
	  }

	  for (int i = 0; i < lockedpoints.Size(); i++)
	  {
		points[lockedpoints[i]].SetType(POINTTYPE.FIXEDPOINT);
	  }


	  /*
	    for (i = 0; i < openelements.Size(); i++)
	    {
	    const Element2d & sel = openelements[i];
	    for (j = 0; j < sel.GetNP(); j++)
	    {
	    INDEX_2 i2;
	    i2.I1() = sel.PNumMod(j+1);
	    i2.I2() = sel.PNumMod(j+2);
	    i2.Sort();
	    boundaryedges->Set (i2, 1);
  
	    points[sel[j]].SetType(FIXEDPOINT);
	    }
	    }
	  */

	  // eltyps.SetSize (GetNE());
	  // eltyps = FREEELEMENT;

	  for (int i = 0; i < GetNSeg(); i++)
	  {
		  Segment seg = segments[i];
		  INDEX_2 i2 = new INDEX_2(seg[0], seg[1]);
		  i2.Sort();

		  //boundaryedges -> Set (i2, 2);
		  segmentht.Set(i2, i);
	  }
	}

	/// additional (temporarily) fix points 
	public void FixPoints(BitArray fixpoints)
	{
	  if (fixpoints.Size() != GetNP())
	  {
		  cerr << "Mesh::FixPoints: sizes don't fit" << "\n";
		  return;
	  }
	  int np = GetNP();
	  for (int i = 1; i <= np; i++)
	  {
		if (fixpoints.Test(i))
		{
			points.Elem(i).SetType(POINTTYPE.FIXEDPOINT);
		}
	  }
	}

	/**
	   finds elements without neighbour and
	   boundary elements without inner element.
	   Results are stored in openelements.
	   if dom == 0, all sub-domains, else subdomain dom */
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
//C++ TO C# CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
	Timer FindOpenElements_t("Mesh::FindOpenElements");

	public void FindOpenElements(int dom = 0)
	{
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	  static Timer t("Mesh::FindOpenElements");
	  RegionTimer reg = new RegionTimer(FindOpenElements_t);

	  int np = GetNP();
	  int ne = GetNE();
	  int nse = GetNSE();

	  Array<int,PointIndex.BASE> numonpoint = new Array<int,PointIndex.BASE>(np);

	  numonpoint = 0;

	  for (ElementIndex ei = 0; ei < ne; ei++)
	  {
		  Element el = this[ei];
		  if (dom == 0 || dom == el.GetIndex())
		  {
			  if (el.GetNP() == 4)
			  {
				  INDEX_4 i4 = new INDEX_4(el[0], el[1], el[2], el[3]);
				  i4.Sort();
				  numonpoint[i4.I1()]++;
				  numonpoint[i4.I2()]++;
			  }
			  else
			  {
				for (int j = 0; j < el.GetNP(); j++)
				{
				  numonpoint[el[j]]++;
				}
			  }
		  }
	  }

	  TABLE<ElementIndex,PointIndex.BASE> elsonpoint = new TABLE<ElementIndex,PointIndex.BASE>(numonpoint);
	  for (ElementIndex ei = 0; ei < ne; ei++)
	  {
		  Element el = this[ei];
		  if (dom == 0 || dom == el.GetIndex())
		  {
			  if (el.GetNP() == 4)
			  {
				  INDEX_4 i4 = new INDEX_4(el[0], el[1], el[2], el[3]);
				  i4.Sort();
				  elsonpoint.Add(i4.I1(), ei);
				  elsonpoint.Add(i4.I2(), ei);
			  }
			  else
			  {
				for (int j = 0; j < el.GetNP(); j++)
				{
				  elsonpoint.Add(el[j], ei);
				}
			  }
		  }
	  }


	  Array<char, 1> hasface = new Array<char, 1>(GetNFD());

	  int i;
	  for (i = 1; i <= GetNFD(); i++)
	  {
		  int domin = GetFaceDescriptor(i).DomainIn();
		  int domout = GetFaceDescriptor(i).DomainOut();
		  hasface[i] = (dom == 0 && (domin != 0 || domout != 0)) || (dom != 0 && (domin == dom || domout == dom));
	  }

	  numonpoint = 0;
	  for (SurfaceElementIndex sii = 0; sii < nse; sii++)
	  {
		  int ind = surfelements[sii].GetIndex();
		  /*
		    if (
		    GetFaceDescriptor(ind).DomainIn() &&
		    (dom == 0 || dom == GetFaceDescriptor(ind).DomainIn())
		    ||
		    GetFaceDescriptor(ind).DomainOut() &&
		    (dom == 0 || dom == GetFaceDescriptor(ind).DomainOut())
		    )
		  */
		  if (hasface[ind])
		  {
			  /*
			    Element2d hel = surfelements[i];
			    hel.NormalizeNumbering();
			    numonpoint[hel[0]]++;
			  */
			  Element2d hel = surfelements[sii];
			  int mini = 0;
			  for (int j = 1; j < hel.GetNP(); j++)
			  {
				if (hel[j] < hel[mini])
				{
				  mini = j;
				}
			  }
			  numonpoint[hel[mini]]++;
		  }
	  }

	  TABLE<SurfaceElementIndex,PointIndex.BASE> selsonpoint = new TABLE<SurfaceElementIndex,PointIndex.BASE>(numonpoint);
	  for (SurfaceElementIndex sii = 0; sii < nse; sii++)
	  {
		  int ind = surfelements[sii].GetIndex();

		  /*
		    if (
		    GetFaceDescriptor(ind).DomainIn() &&
		    (dom == 0 || dom == GetFaceDescriptor(ind).DomainIn())
		    ||
		    GetFaceDescriptor(ind).DomainOut() &&
		    (dom == 0 || dom == GetFaceDescriptor(ind).DomainOut())
		    )
		  */
		  if (hasface[ind])
		  {
			  /*
			    Element2d hel = surfelements[i];
			    hel.NormalizeNumbering();
			    selsonpoint.Add (hel[0], i);
			  */
			  Element2d hel = surfelements[sii];
			  int mini = 0;
			  for (int j = 1; j < hel.GetNP(); j++)
			  {
				if (hel[j] < hel[mini])
				{
				  mini = j;
				}
			  }
			  selsonpoint.Add(hel[mini], sii);
		  }
	  }


	  int ii;
	  PointIndex pi = new PointIndex();
	  SurfaceElementIndex sei = new SurfaceElementIndex();
	  // Element2d hel;


	  INDEX_3_CLOSED_HASHTABLE<INDEX_2> faceht = new INDEX_3_CLOSED_HASHTABLE<INDEX_2>(100);
	  openelements.SetSize(0);

	  for (PointIndex pi = points.Begin(); pi < points.End(); pi++)
	  {
		if ((selsonpoint[pi].Size() + elsonpoint[pi].Size()) != 0)
		{
			faceht.SetSize(2 * selsonpoint[pi].Size() + 4 * elsonpoint[pi].Size());

			FlatArray<SurfaceElementIndex> row = selsonpoint[pi];
			for (ii = 0; ii < row.Size(); ii++)
			{
				Element2d hel = SurfaceElement(row[ii]);
				if (hel.GetType() == ELEMENT_TYPE.TRIG6)
				{
					hel.SetType(ELEMENT_TYPE.TRIG);
				}
				int ind = hel.GetIndex();

				if ((GetFaceDescriptor(ind).DomainIn()) != 0 && (dom == 0 || dom == GetFaceDescriptor(ind).DomainIn()))
				{
					hel.NormalizeNumbering();
					if (hel.PNum(1) == pi)
					{
						INDEX_3 i3 = new INDEX_3(hel[0], hel[1], hel[2]);
						INDEX_2 i2(GetFaceDescriptor(ind).DomainIn(), (hel.GetNP() == 3) ? new PointIndex(PointIndex.BASE-1) : hel.PNum(4));
						faceht.Set(i3, i2);
					}
				}
				if ((GetFaceDescriptor(ind).DomainOut()) != 0 && (dom == 0 || dom == GetFaceDescriptor(ind).DomainOut()))
				{
					hel.Invert();
					hel.NormalizeNumbering();
					if (hel.PNum(1) == pi)
					{
						INDEX_3 i3 = new INDEX_3(hel[0], hel[1], hel[2]);
						INDEX_2 i2(GetFaceDescriptor(ind).DomainOut(), (hel.GetNP() == 3) ? new PointIndex(PointIndex.BASE-1) : hel.PNum(4));
						faceht.Set(i3, i2);
					}
				}
			}


			FlatArray<ElementIndex> rowel = elsonpoint[pi];
			for (ii = 0; ii < rowel.Size(); ii++)
			{
				Element el = VolumeElement(rowel[ii]);

				if (dom == 0 || el.GetIndex() == dom)
				{
					for (int j = 1; j <= el.GetNFaces(); j++)
					{
						Element2d hel = new Element2d(ELEMENT_TYPE.TRIG);
						el.GetFace(j, hel);
						hel.Invert();
						hel.NormalizeNumbering();

						if (hel[0] == pi)
						{
							INDEX_3 i3 = new INDEX_3(hel[0], hel[1], hel[2]);

							if (faceht.Used(i3))
							{
								INDEX_2 i2 = faceht.Get(i3);
								if (i2.I1() == el.GetIndex())
								{
									i2.I1() = PointIndex.BASE-1;
									faceht.Set(i3, i2);
								}
								else
								{
									if (i2.I1() == 0)
									{
										PrintSysError("more elements on face");
										(*testout) << "more elements on face!!!" << "\n";
										(*testout) << "el = " << el << "\n";
										(*testout) << "hel = " << hel << "\n";
										(*testout) << "face = " << i3 << "\n";
										(*testout) << "points = " << "\n";
										for (int jj = 1; jj <= 3; jj++)
										{
										  (*testout) << "p = " << new Point(i3.I(jj)) << "\n";
										}
									}
								}
							}
							else
							{
								hel.Invert();
								hel.NormalizeNumbering();
								INDEX_3 i3 = new INDEX_3(hel[0], hel[1], hel[2]);
								INDEX_2 i2(el.GetIndex(), (hel.GetNP() == 3) ? new PointIndex(PointIndex.BASE-1) : hel[3]);
								faceht.Set(i3, i2);
							}
						}
					}
				}
			}
			for (int i = 0; i < faceht.Size(); i++)
			{
			  if (faceht.UsedPos(i))
			  {
				  INDEX_3 i3 = new INDEX_3();
				  INDEX_2 i2 = new INDEX_2();
				  faceht.GetData(i, ref i3, ref i2);
				  if (i2.I1() != PointIndex.BASE-1)
				  {
					  // Element2d tri;
					  // tri.SetType ( (i2.I2() == PointIndex::BASE-1) ? TRIG : QUAD);
					  Element2d tri((i2.I2() == PointIndex.BASE-1) ? ELEMENT_TYPE.TRIG : ELEMENT_TYPE.QUAD);
					  for (int l = 0; l < 3; l++)
					  {
						tri[l] = i3.I(l + 1);
					  }
					  tri.PNum(4) = i2.I2();
					  tri.SetIndex(i2.I1());

					  //	tri.Invert();

					  openelements.Append(tri);
				  }
			  }
			}
		}
	  }

	  int cnt3 = 0;
	  for (i = 0; i < openelements.Size(); i++)
	  {
		if (openelements[i].GetNP() == 3)
		{
		  cnt3++;
		}
	  }

	  int cnt4 = openelements.Size() - cnt3;


	  MyStr treequad = new MyStr();
	  if (cnt4 != 0)
	  {
		treequad = new MyStr(" (") + new MyStr(cnt3) + new MyStr(" + ") + new MyStr(cnt4) + new MyStr(")");
	  }

	  PrintMessage(5, openelements.Size(), treequad.functorMethod, " open elements");

	  BuildBoundaryEdges();


	  for (int i = 1; i <= openelements.Size(); i++)
	  {
		  Element2d sel = openelements.Get(i);

		  if (boundaryedges != null)
		  {
			for (int j = 1; j <= sel.GetNP(); j++)
			{
				INDEX_2 i2 = new INDEX_2();
				i2.I1() = sel.PNumMod(j);
				i2.I2() = sel.PNumMod(j + 1);
				i2.Sort();
				boundaryedges.Set(i2, 1);
			}
		  }

		  for (int j = 1; j <= 3; j++)
		  {
			  PointIndex pi = sel.PNum(j);
			  if (pi < points.End())
			  {
				points[pi].SetType(POINTTYPE.FIXEDPOINT);
			  }
		  }
	  }



	  /*
	    for (i = 1; i <= GetNSeg(); i++)
	    {
	    const Segment & seg = LineSegment(i);
	    INDEX_2 i2(seg[0], seg[1]);
	    i2.Sort();
  
	    if (!boundaryedges->Used (i2))
	    cerr << "WARNING: no boundedge, but seg edge: " << i2 << endl;
  
	    boundaryedges -> Set (i2, 2);
	    segmentht -> Set (i2, i-1);
	    }
	  */
	}


	/**
	   finds segments without surface element,
	   and surface elements without neighbours.
	   store in opensegmentsy
	*/
	public void FindOpenSegments(int surfnr = 0)
	{
	  // int i, j, k;

	  // new version, general elements
	  // hash index: pnum1-2
	  // hash data : surfnr,  surfel-nr (pos) or segment nr(neg)
	  INDEX_2_HASHTABLE<INDEX_2> faceht = new INDEX_2_HASHTABLE<INDEX_2>(4 * GetNSE() + GetNSeg() + 1);

	  PrintMessage(5, "Test Opensegments");
	  for (int i = 1; i <= GetNSeg(); i++)
	  {
		  Segment seg = LineSegment(i);

		  if (surfnr == 0 || seg.si == surfnr)
		  {
			  INDEX_2 key = new INDEX_2(seg[0], seg[1]);
			  INDEX_2 data = new INDEX_2(seg.si, -i);

			  if (faceht.Used(key))
			  {
				  cerr << "ERROR: Segment " << seg << " already used" << "\n";
				  (*testout) << "ERROR: Segment " << seg << " already used" << "\n";
			  }

			  faceht.Set(key, data);
		  }
	  }


	  for (int i = 1; i <= GetNSeg(); i++)
	  {
		  Segment seg = LineSegment(i);

		  if (surfnr == 0 || seg.si == surfnr)
		  {
			  INDEX_2 key = new INDEX_2(seg[1], seg[0]);
			  if (!faceht.Used(key))
			  {
				  cerr << "ERROR: Segment " << seg << " brother not used" << "\n";
				  (*testout) << "ERROR: Segment " << seg << " brother not used" << "\n";
			  }
		  }
	  }

	  // bool buggy = false;
	  // ofstream bout("buggy.out");

	  for (int i = 1; i <= GetNSE(); i++)
	  {
		  Element2d el = SurfaceElement(i);
		  if (el.IsDeleted())
		  {
			  continue;
		  }

		  if (surfnr == 0 || el.GetIndex() == surfnr)
		  {
			  for (int j = 1; j <= el.GetNP(); j++)
			  {
				  INDEX_2 seg = new INDEX_2(el.PNumMod(j), el.PNumMod(j + 1));
				  INDEX_2 data = new INDEX_2();

				  if (seg.I1() < PointIndex.BASE || seg.I2() < PointIndex.BASE)
				  {
					cerr << "seg = " << seg << "\n";
				  }

				  if (faceht.Used(seg))
				  {
					  data = faceht.Get(seg);
					  if (data.I1() == el.GetIndex())
					  {
						  data.I1() = 0;
						  faceht.Set(seg, data);
					  }
					  else
					  {
			  // buggy = true;
						  PrintWarning("hash table si not fitting for segment: ", seg.I1(), "-", seg.I2(), " other = ", data.I2());
			  // cout << "me: index = " << el.GetIndex() << ", el = " << el << endl;

			  /*
			  bout << "has index = " << seg << endl;
			  bout << "hash value = " << faceht.HashValue (seg) << endl;
  
			  if (data.I2() > 0)
			    {
			      int io = data.I2();
			      cout << "other trig: index = " << SurfaceElement(io).GetIndex()
			  	 << ", el = " << SurfaceElement(io) << endl;
			    }
			  else
			    {
			      cout << "other seg " << -data.I2() << ", si = " << data.I1() << endl;
			    }
  
			  
			  bout << "me: index = " << el.GetIndex() << ", el = " << el << endl;
			  if (data.I2() > 0)
			    {
			      int io = data.I2();
			      bout << "other trig: index = " << SurfaceElement(io).GetIndex()
			  	 << ", el = " << SurfaceElement(io) << endl;
			    }
			  else
			    {
			      bout << "other seg " << -data.I2() << ", si = " << data.I1() << endl;
			    }
			  */
					  }
				  }
				  else
				  {
					  netgen.GlobalMembers.Swap(ref seg.I1(), ref seg.I2());
					  data.I1() = el.GetIndex();
					  data.I2() = i;

					  faceht.Set(seg, data);
				  }
			  }
		  }
	  }

	  /*
	  if (buggy)
	    {
	  for (int i = 1; i <= GetNSeg(); i++)
		bout << "seg" << i << " " << LineSegment(i) << endl;
  
	  for (int i = 1; i <= GetNSE(); i++)
		bout << "sel" << i << " " << SurfaceElement(i) << " ind = "
			 << SurfaceElement(i).GetIndex() << endl;
  
	  bout << "hashtable: " << endl;
	  for (int j = 1; j <= faceht.GetNBags(); j++)
		{
		  bout << "bag " << j << ":" << endl;
		  for (int k = 1; k <= faceht.GetBagSize(j); k++)
			{
		  INDEX_2 i2, data;
		  faceht.GetData (j, k, i2, data);
		  bout << "key = " << i2 << ", data = " << data << endl;
			}
		}
	  exit(1);
	    }
	  */

	  (*testout) << "open segments: " << "\n";
	  opensegments.SetSize(0);
	  for (int i = 1; i <= faceht.GetNBags(); i++)
	  {
		for (int j = 1; j <= faceht.GetBagSize(i); j++)
		{
			INDEX_2 i2 = new INDEX_2();
			INDEX_2 data = new INDEX_2();
			faceht.GetData(i, j, ref i2, ref data);
			if (data.I1() != null) // surfnr
			{
				Segment seg = new Segment();
				seg[0] = i2.I1();
				seg[1] = i2.I2();
				seg.si = data.I1();

				// find geomdata:
				if (data.I2() > 0)
				{
					// segment due to triangle
					Element2d el = SurfaceElement(data.I2());
					for (int k = 1; k <= el.GetNP(); k++)
					{
						if (seg[0] == el.PNum(k))
						{
						  seg.geominfo[0] = el.GeomInfoPi(k);
						}
						if (seg[1] == el.PNum(k))
						{
						  seg.geominfo[1] = el.GeomInfoPi(k);
						}
					}

					(*testout) << "trig seg: ";
				}
				else
				{
					// segment due to line
					Segment lseg = LineSegment(-data.I2());
					seg.geominfo[0] = lseg.geominfo[0];
					seg.geominfo[1] = lseg.geominfo[1];

					(*testout) << "line seg: ";
				}

				(*testout) << seg[0] << " - " << seg[1] << " len = " << netgen.GlobalMembers.Dist(new Point(seg[0]), new Point(seg[1])) << "\n";

				opensegments.Append(seg);
				if (seg.geominfo[0].trignum <= 0 || seg.geominfo[1].trignum <= 0)
				{
					(*testout) << "Problem with open segment: " << seg << "\n";
				}

			}
		}
	  }

	  PrintMessage(3, opensegments.Size(), " open segments found");
	  (*testout) << opensegments.Size() << " open segments found" << "\n";

	  /*
	    ptyps.SetSize (GetNP());
	    for (i = 1; i <= ptyps.Size(); i++)
	    ptyps.Elem(i) = SURFACEPOINT;
  
	    for (i = 1; i <= GetNSeg(); i++)
	    {
	    const Segment & seg = LineSegment (i);
	    ptyps.Elem(seg[0]) = EDGEPOINT;
	    ptyps.Elem(seg[1]) = EDGEPOINT;
	    }
	    for (i = 1; i <= GetNOpenSegments(); i++)
	    {
	    const Segment & seg = GetOpenSegment (i);
	    ptyps.Elem(seg[0]) = EDGEPOINT;
	    ptyps.Elem(seg[1]) = EDGEPOINT;
	    }
	  */
	  for (int i = 1; i <= points.Size(); i++)
	  {
		points.Elem(i).SetType(POINTTYPE.SURFACEPOINT);
	  }

	  for (int i = 1; i <= GetNSeg(); i++)
	  {
		  Segment seg = LineSegment(i);
		  points[seg[0]].SetType(POINTTYPE.EDGEPOINT);
		  points[seg[1]].SetType(POINTTYPE.EDGEPOINT);
	  }
	  for (int i = 1; i <= GetNOpenSegments(); i++)
	  {
		  Segment seg = GetOpenSegment(i);
		  points[seg[0]].SetType(POINTTYPE.EDGEPOINT);
		  points[seg[1]].SetType(POINTTYPE.EDGEPOINT);
	  }



	  /*
  
	  for (i = 1; i <= openelements.Size(); i++)
	  {
	  const Element2d & sel = openelements.Get(i);
  
	  if (boundaryedges)
	  for (j = 1; j <= sel.GetNP(); j++)
	  {
	  INDEX_2 i2;
	  i2.I1() = sel.PNumMod(j);
	  i2.I2() = sel.PNumMod(j+1);
	  i2.Sort();
	  boundaryedges->Set (i2, 1);
	  }
  
	  for (j = 1; j <= 3; j++)
	  {
	  int pi = sel.PNum(j);
	  if (pi <= ptyps.Size())
	  ptyps.Elem(pi) = FIXEDPOINT;
	  }
	  }
	  */
	}

	/**
	   remove one layer of surface elements
	*/
	public void RemoveOneLayerSurfaceElements()
	{
	  int np = GetNP();

	  FindOpenSegments();
	  BitArray frontpoints = new BitArray(np + 1); // for 0- and 1-based
	  frontpoints.Clear();

	  for (int i = 1; i <= GetNOpenSegments(); i++)
	  {
		  Segment seg = GetOpenSegment(i);
		  frontpoints.Set(new netgen.Segment(seg[0]));
		  frontpoints.Set(new netgen.Segment(seg[1]));
	  }

	  for (int i = 1; i <= GetNSE(); i++)
	  {
		  Element2d sel = surfelements.Elem(i);
		  bool remove = false;
		  for (int j = 1; j <= sel.GetNP(); j++)
		  {
			if (frontpoints.Test(sel.PNum(j)))
			{
			  remove = true;
			}
		  }
		  if (remove)
		  {
			sel.PNum(1).Invalidate();
		  }
	  }

	  for (int i = surfelements.Size(); i >= 1; i--)
	  {
		  if (!surfelements.Elem(i).PNum(1).IsValid())
		  {
			  surfelements.Elem(i) = surfelements.Last();
			  surfelements.DeleteLast();
		  }
	  }

	  RebuildSurfaceElementLists();
	  /*
	  for (int i = 0; i < facedecoding.Size(); i++)
	    facedecoding[i].firstelement = -1;
	  for (int i = surfelements.Size()-1; i >= 0; i--)
	    {
	      int ind = surfelements[i].GetIndex();
	      surfelements[i].next = facedecoding[ind-1].firstelement;
	      facedecoding[ind-1].firstelement = i;
	    }
	  */

	  timestamp = netgen.GlobalMembers.NextTimeStamp();
	  //  Compress();
	}


	public int GetNOpenSegments()
	{
		return opensegments.Size();
	}
	public Segment GetOpenSegment(int nr)
	{
		return opensegments.Get(nr);
	}

	/**
	   Checks overlap of boundary
	   return == 1, iff overlap
	*/
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
//C++ TO C# CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
	Timer CheckOverlappingBoundary_t("Mesh::CheckOverlappingBoundary");

	public int CheckOverlappingBoundary()
	{
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	  static Timer t("Mesh::CheckOverlappingBoundary");
	  RegionTimer reg = new RegionTimer(CheckOverlappingBoundary_t);

	  int i;
	  int j;
	  int k;

	  Point3d pmin = new Point3d();
	  Point3d pmax = new Point3d();
	  GetBox(ref pmin, ref pmax);
	  BoxTree < 3> setree(pmin, pmax);
	  Array<int> inters = new Array<int>();

	  bool overlap = false;
	  bool incons_layers = false;


	  for (i = 1; i <= GetNSE(); i++)
	  {
		SurfaceElement(i).badel = 0;
	  }


	  for (i = 1; i <= GetNSE(); i++)
	  {
		  Element2d tri = SurfaceElement(i);

		  Point3d tpmin = new Point3d(new Point(tri[0]));
		  Point3d tpmax = new Point3d(tpmin);

		  for (k = 1; k < tri.GetNP(); k++)
		  {
			  tpmin.SetToMin(new Point(tri[k]));
			  tpmax.SetToMax(new Point(tri[k]));
		  }
		  Vec3d diag = new Vec3d(tpmin, tpmax);

		  tpmax.CopyFrom(tpmax + 0.1 * diag);
		  tpmin.CopyFrom(tpmin - 0.1 * diag);

		  setree.Insert(tpmin, tpmax, i);
	  }

	  for (i = 1; i <= GetNSE(); i++)
	  {
		  Element2d tri = SurfaceElement(i);

		  Point3d tpmin = new Point3d(new Point(tri[0]));
		  Point3d tpmax = new Point3d(tpmin);

		  for (k = 1; k < tri.GetNP(); k++)
		  {
			  tpmin.SetToMin(new Point(tri[k]));
			  tpmax.SetToMax(new Point(tri[k]));
		  }

		  setree.GetIntersecting(tpmin, tpmax, inters);

		  for (j = 1; j <= inters.Size(); j++)
		  {
			  Element2d tri2 = SurfaceElement(inters.Get(j));

			  if (this[tri[0]].GetLayer() != this[tri2[0]].GetLayer())
			  {
				continue;
			  }

			  if (this[tri[0]].GetLayer() != this[tri[1]].GetLayer() || this[tri[0]].GetLayer() != this[tri[2]].GetLayer())
			  {
				  incons_layers = true;
				  Console.Write("inconsistent layers in triangle");
				  Console.Write("\n");
			  }


			  const netgen.Point < 3> *trip1[3], *trip2[3];
			  for (k = 1; k <= 3; k++)
			  {
				  trip1[k - 1] = &new Point(tri.PNum(k));
				  trip2[k - 1] = &new Point(tri2.PNum(k));
			  }

			  if (IntersectTriangleTriangle(trip1[0], trip2[0]))
			  {
				  overlap = true;
				  PrintWarning("Intersecting elements ",i, " and ", inters.Get(j));

				  (*testout) << "Intersecting: " << "\n";
				  (*testout) << "openelement " << i << " with open element " << inters.Get(j) << "\n";

				  Console.Write("el1 = ");
				  Console.Write(tri);
				  Console.Write("\n");
				  Console.Write("el2 = ");
				  Console.Write(tri2);
				  Console.Write("\n");
				  Console.Write("layer1 = ");
				  Console.Write(this[tri[0]].GetLayer());
				  Console.Write("\n");
				  Console.Write("layer2 = ");
				  Console.Write(this[tri2[0]].GetLayer());
				  Console.Write("\n");


				  for (k = 1; k <= 3; k++)
				  {
					(*testout) << tri.PNum(k) << "  ";
				  }
				  (*testout) << "\n";
				  for (k = 1; k <= 3; k++)
				  {
					(*testout) << tri2.PNum(k) << "  ";
				  }
				  (*testout) << "\n";

				  for (k = 0; k <= 2; k++)
				  {
					(*testout) << *trip1[k] << "   ";
				  }
				  (*testout) << "\n";
				  for (k = 0; k <= 2; k++)
				  {
					(*testout) << *trip2[k] << "   ";
				  }
				  (*testout) << "\n";

				  (*testout) << "Face1 = " << GetFaceDescriptor(tri.GetIndex()) << "\n";
				  (*testout) << "Face1 = " << GetFaceDescriptor(tri2.GetIndex()) << "\n";

				  /*
				    INDEX_3 i3(tri.PNum(1), tri.PNum(2), tri.PNum(3));
				    i3.Sort();
				    for (k = 1; k <= GetNSE(); k++)
				    {
				    const Element2d & el2 = SurfaceElement(k);
				    INDEX_3 i3b(el2.PNum(1), el2.PNum(2), el2.PNum(3));
				    i3b.Sort();
				    if (i3 == i3b)
				    {
				    SurfaceElement(k).badel = 1;
				    }
				    }
				  */
				  SurfaceElement(i).badel = 1;
				  SurfaceElement(inters.Get(j)).badel = 1;
			  }
		  }
	  }

	  // bug 'fix'
	  if (incons_layers)
	  {
		  overlap = false;
	  }

	  return overlap;
	}

	/**
	   Checks consistent boundary
	   return == 0, everything ok
	*/
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int CheckConsistentBoundary() const
	public int CheckConsistentBoundary()
	{
	  int nf = GetNOpenElements();
	  INDEX_2_HASHTABLE<int> edges = new INDEX_2_HASHTABLE<int>(nf + 2);
	  INDEX_2 i2 = new INDEX_2();
	  INDEX_2 i2s = new INDEX_2();
	  INDEX_2 edge = new INDEX_2();
	  int err = 0;

	  for (int i = 1; i <= nf; i++)
	  {
		  Element2d sel = OpenElement(i);

		  for (int j = 1; j <= sel.GetNP(); j++)
		  {
			  i2.I1() = sel.PNumMod(j);
			  i2.I2() = sel.PNumMod(j + 1);

			  int sign = (i2.I2() > i2.I1()) ? 1 : -1;
			  i2.Sort();
			  if (!edges.Used(i2))
			  {
				edges.Set(i2, 0);
			  }
			  edges.Set(i2, edges.Get(i2) + sign);
		  }
	  }

	  for (int i = 1; i <= edges.GetNBags(); i++)
	  {
		for (int j = 1; j <= edges.GetBagSize(i); j++)
		{
			int cnt = 0;
			edges.GetData(i, j, ref i2, ref cnt);
			if (cnt != 0)
			{
				PrintError("Edge ", i2.I1(), " - ", i2.I2(), " multiple times in surface mesh");

				(*testout) << "Edge " << i2 << " multiple times in surface mesh" << "\n";
//C++ TO C# CONVERTER TODO TASK: The following line was determined to be a copy assignment (rather than a reference assignment) - this should be verified and a 'CopyFrom' method should be created:
//ORIGINAL LINE: i2s = i2;
				i2s.CopyFrom(i2);
				i2s.Sort();
				for (int k = 1; k <= nf; k++)
				{
					Element2d sel = OpenElement(k);
					for (int l = 1; l <= sel.GetNP(); l++)
					{
						edge.I1() = sel.PNumMod(l);
						edge.I2() = sel.PNumMod(l + 1);
						edge.Sort();

						if (edge == i2s)
						{
						  (*testout) << "edge of element " << sel << "\n";
						}
					}
				}


				err = 2;
			}
		}
	  }

	  return err;
	}

	/*
	  checks element orientation
	*/
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int CheckVolumeMesh() const
	public int CheckVolumeMesh()
	{
	  PrintMessage(3, "Checking volume mesh");

	  int ne = GetNE();
	  DenseMatrix dtrans = new DenseMatrix(3, 3);
	  int i;
	  int j;

	  PrintMessage(5, "elements: ", ne);
	  for (i = 1; i <= ne; i++)
	  {
		  Element el = (Element) VolumeElement(i);
		  el.flags.badel = 0;
		  int nip = el.GetNIP();
		  for (j = 1; j <= nip; j++)
		  {
			  el.GetTransformation(j, Points(), dtrans.functorMethod);
			  double det = dtrans.Det();
			  if (det > 0)
			  {
				  PrintError("Element ", i, " has wrong orientation");
				  el.flags.badel = 1;
			  }
		  }
	  }

	  return 0;
	}


	/**
	   finds average h of surface surfnr if surfnr > 0,
	   else of all surfaces.
	*/
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: double AverageH(int surfnr = 0) const
	public double AverageH(int surfnr = 0)
	{
	  int i;
	  int j;
	  int n;
	  double hi;
	  double hsum;
	  double maxh = 0;
	  double minh = 1e10;

	  hsum = 0;
	  n = 0;
	  for (i = 1; i <= GetNSE(); i++)
	  {
		  Element2d el = SurfaceElement(i);
		  if (surfnr == 0 || el.GetIndex() == surfnr)
		  {
			  for (j = 1; j <= 3; j++)
			  {
				  hi = netgen.GlobalMembers.Dist(new Point(el.PNumMod(j)), new Point(el.PNumMod(j + 1)));

				  hsum += hi;

				  if (hi > maxh)
				  {
					  maxh = hi;
				  }
				  if (hi < minh)
				  {
					  minh = hi;
				  }
				  n++;
			  }
		  }
	  }

	  PrintMessage(5, "minh = ", minh, " avh = ", (hsum / n), " maxh = ", maxh);
	  return (hsum / n);
	}

	/// Calculates localh 
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
//C++ TO C# CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
	Timer CalcLocalH_t("Mesh::CalcLocalH");

	public void CalcLocalH(double grading)
	{
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	  static Timer t("Mesh::CalcLocalH");
	  RegionTimer reg = new RegionTimer(CalcLocalH_t);

	  if (lochfunc == null)
	  {
		  Point3d pmin = new Point3d();
		  Point3d pmax = new Point3d();
		  GetBox(ref pmin, ref pmax);
		  // SetLocalH (pmin, pmax, mparam.grading);
	  SetLocalH(new netgen.Point3d(pmin), new netgen.Point3d(pmax), grading);
	  }

	  PrintMessage(3, "CalcLocalH: ", GetNP(), " Points ", GetNE(), " Elements ", GetNSE(), " Surface Elements");


	  for (int i = 0; i < GetNSE(); i++)
	  {
		  Element2d el = surfelements[i];
		  int j;

		  if (el.GetNP() == 3)
		  {
			  double hel = -1;
			  for (j = 1; j <= 3; j++)
			  {
				  Point3d p1 = points[el.PNumMod(j)];
				  Point3d p2 = points[el.PNumMod(j + 1)];

				  /*
				    INDEX_2 i21(el.PNumMod(j), el.PNumMod(j+1));
				    INDEX_2 i22(el.PNumMod(j+1), el.PNumMod(j));
				    if (! identifiedpoints->Used (i21) &&
				    ! identifiedpoints->Used (i22) )
				  */
				  if (!ident.UsedSymmetric(el.PNumMod(j), el.PNumMod(j + 1)))
				  {
					  double hedge = netgen.GlobalMembers.Dist(p1, p2);
					  if (hedge > hel)
					  {
						hel = hedge;
					  }
					  //		  lochfunc->SetH (Center (p1, p2), 2 * Dist (p1, p2));
					  //		  (*testout) << "trigseth, p1,2 = " << el.PNumMod(j) << ", " << el.PNumMod(j+1)
					  //			     << " h = " << (2 * Dist(p1, p2)) << endl;
				  }
			  }

			  if (hel > 0)
			  {
				  Point3d p1 = points[el.PNum(1)];
				  Point3d p2 = points[el.PNum(2)];
				  Point3d p3 = points[el.PNum(3)];
				  lochfunc.SetH(netgen.GlobalMembers.Center(p1, p2, p3), hel);
			  }
		  }
		  else
		  {
			  {
				Point3d p1 = points[el.PNum(1)];
				Point3d p2 = points[el.PNum(2)];
				lochfunc.SetH(netgen.GlobalMembers.Center(p1, p2), 2 * netgen.GlobalMembers.Dist(p1, p2));
			  }
			  {
				Point3d p1 = points[el.PNum(3)];
				Point3d p2 = points[el.PNum(4)];
				lochfunc.SetH(netgen.GlobalMembers.Center(p1, p2), 2 * netgen.GlobalMembers.Dist(p1, p2));
			  }
		  }
	  }

	  for (int i = 0; i < GetNSeg(); i++)
	  {
		  Segment seg = segments[i];
		  Point3d p1 = points[seg[0]];
		  Point3d p2 = points[seg[1]];
		  /*
		    INDEX_2 i21(seg[0], seg[1]);
		    INDEX_2 i22(seg[1], seg[0]);
		    if (identifiedpoints)
		    if (!identifiedpoints->Used (i21) && !identifiedpoints->Used (i22))
		  */
		  if (!ident.UsedSymmetric(new netgen.Segment(seg[0]), new netgen.Segment(seg[1])))
		  {
			  lochfunc.SetH(netgen.GlobalMembers.Center(p1, p2), netgen.GlobalMembers.Dist(p1, p2));
		  }
	  }
	  /*
	    cerr << "do vol" << endl;
	    for (i = 1; i <= GetNE(); i++)
	    {
	    const Element & el = VolumeElement(i);
	    if (el.GetType() == TET)
	    {
	    int j, k;
	    for (j = 2; j <= 4; j++)
	    for (k = 1; k < j; k++)
	    {
	    const Point3d & p1 = Point (el.PNum(j));
	    const Point3d & p2 = Point (el.PNum(k));
	    lochfunc->SetH (Center (p1, p2), 2 * Dist (p1, p2));
	    (*testout) << "set vol h to " << (2 * Dist (p1, p2)) << endl;
  
	    }
	    }
	    }
	  */

	  /*
	    const char * meshsizefilename =
	    globflags.GetStringFlag ("meshsize", NULL);
	    if (meshsizefilename)
	    {
	    ifstream msf(meshsizefilename);
	    if (msf)
	    {
	    int nmsp;
	    msf >> nmsp;
	    for (i = 1; i <= nmsp; i++)
	    {
	    Point3d pi;
	    double hi;
	    msf >> pi.X() >> pi.Y() >> pi.Z();
	    msf >> hi;
	    lochfunc->SetH (pi, hi);
	    }
	    }
	    }
	  */
	  //  lochfunc -> Convexify();
	  //  lochfunc -> PrintMemInfo (cout);
	}

	///
	public void SetLocalH(netgen.Point < 3> pmin, netgen.Point < 3> pmax, double grading)
	{
	  Point < 3> c = netgen.GlobalMembers.Center(pmin.functorMethod, pmax.functorMethod);
	  double d = netgen.GlobalMembers.max3(pmax.functorMethod(0) - pmin.functorMethod(0), pmax.functorMethod(1) - pmin.functorMethod(1), pmax.functorMethod(2) - pmin.functorMethod(2));
	  d /= 2;
	  Point < 3> pmin2 = c - Vec < 3> (d, d, d);
	  Point < 3> pmax2 = c + Vec < 3> (d, d, d);

	  if (lochfunc != null)
	  {
		  lochfunc.Dispose();
	  }
	  lochfunc = new LocalH(pmin2, pmax2, grading, dimension);
	}

	///
	public void RestrictLocalH(Point3d p, double hloc)
	{
	  if (hloc < hmin)
	  {
		hloc = hmin;
	  }

	  //cout << "restrict h in " << p << " to " << hloc << endl;
	  if (lochfunc == null)
	  {
		  PrintWarning("RestrictLocalH called, creating mesh-size tree");

		  Point3d boxmin = new Point3d();
		  Point3d boxmax = new Point3d();
		  GetBox(ref boxmin, ref boxmax);
		  SetLocalH(new netgen.Point3d(boxmin), new netgen.Point3d(boxmax), 0.8);
	  }

	  lochfunc.SetH(new netgen.Point3d(p), hloc);
	}

	///
	public void RestrictLocalHLine(Point3d p1, Point3d p2, double hloc)
	{
	  if (hloc < hmin)
	  {
		hloc = hmin;
	  }

	  // cout << "restrict h along " << p1 << " - " << p2 << " to " << hloc << endl;
	  int i;
	  int steps = (int)(netgen.GlobalMembers.Dist(p1, p2) / hloc) + 2;
	  Vec3d v = new Vec3d(p1, p2);

	  for (i = 0; i <= steps; i++)
	  {
		  Point3d p = p1 + ((double)i / (double)steps * v);
		  RestrictLocalH(p, hloc);
	  }
	}

	/// number of elements per radius
	public void CalcLocalHFromSurfaceCurvature(double grading, double elperr)
	{
	  PrintMessage(3, "Calculating local h from surface curvature");

	  if (lochfunc == null)
	  {
		  Point3d pmin = new Point3d();
		  Point3d pmax = new Point3d();
		  GetBox(ref pmin, ref pmax);

		  // SetLocalH (pmin, pmax, mparam.grading);
	  SetLocalH(new netgen.Point3d(pmin), new netgen.Point3d(pmax), grading);
	  }


	  INDEX_2_HASHTABLE<int> edges = new INDEX_2_HASHTABLE<int>(3 * GetNP() + 2);
	  INDEX_2_HASHTABLE<int> bedges = new INDEX_2_HASHTABLE<int>(GetNSeg() + 2);
	  int i;
	  int j;

	  for (i = 1; i <= GetNSeg(); i++)
	  {
		  Segment seg = LineSegment(i);
		  INDEX_2 i2 = new INDEX_2(seg[0], seg[1]);
		  i2.Sort();
		  bedges.Set(i2, 1);
	  }
	  for (i = 1; i <= GetNSE(); i++)
	  {
		  Element2d sel = SurfaceElement(i);
		  if (sel.PNum(1) == null)
		  {
			continue;
		  }
		  for (j = 1; j <= 3; j++)
		  {
			  INDEX_2 i2 = new INDEX_2(sel.PNumMod(j), sel.PNumMod(j + 1));
			  i2.Sort();
			  if (bedges.Used(i2))
			  {
				  continue;
			  }

			  if (edges.Used(i2))
			  {
				  int other = edges.Get(i2);

				  Element2d elother = SurfaceElement(other);

				  int pi3 = 1;
				  while ((sel.PNum(pi3) == i2.I1()) || (sel.PNum(pi3) == i2.I2()))
				  {
					pi3++;
				  }
				  pi3 = sel.PNum(pi3);

				  int pi4 = 1;
				  while ((elother.PNum(pi4) == i2.I1()) || (elother.PNum(pi4) == i2.I2()))
				  {
					pi4++;
				  }
				  pi4 = elother.PNum(pi4);

				  double rad = ComputeCylinderRadius(new Point(i2.I1()), new Point(i2.I2()), new Point(pi3), new Point(pi4));

				  RestrictLocalHLine(new Point(i2.I1()), new Point(i2.I2()), rad / elperr);


				  /*
				    (*testout) << "pi1,2, 3, 4 = " << i2.I1() << ", " << i2.I2() << ", " << pi3 << ", " << pi4
				    << " p1 = " << Point(i2.I1())
				    << ", p2 = " << Point(i2.I2())
				    //			 << ", p3 = " << Point(pi3)
				    //			 << ", p4 = " << Point(pi4)
				    << ", rad = " << rad << endl;
				  */
			  }
			  else
			  {
				edges.Set(i2, i);
			  }
		  }
	  }


	  // Restrict h due to line segments

	  for (i = 1; i <= GetNSeg(); i++)
	  {
		  Segment seg = LineSegment(i);
		  Point3d p1 = new Point(seg[0]);
		  Point3d p2 = new Point(seg[1]);
		  RestrictLocalH(netgen.GlobalMembers.Center(p1, p2), netgen.GlobalMembers.Dist(p1, p2));
	  }



	  /*
  
  
	  int i, j;
	  int np = GetNP();
	  int nseg = GetNSeg();
	  int nse = GetNSE();
  
	  Array<Vec3d> normals(np);
	  BitArray linepoint(np);
  
	  linepoint.Clear();
	  for (i = 1; i <= nseg; i++)
	  {
	  linepoint.Set (LineSegment(i)[0]);
	  linepoint.Set (LineSegment(i)[1]);
	  }
  
	  for (i = 1; i <= np; i++)
	  normals.Elem(i) = Vec3d(0,0,0);
  
	  for (i = 1; i <= nse; i++)
	  {
	  Element2d & el = SurfaceElement(i);
	  Vec3d nf = Cross (Vec3d (Point (el.PNum(1)), Point(el.PNum(2))),
	  Vec3d (Point (el.PNum(1)), Point(el.PNum(3))));
	  for (j = 1; j <= 3; j++)
	  normals.Elem(el.PNum(j)) += nf;
	  }
  
	  for (i = 1; i <= np; i++)
	  normals.Elem(i) /= (1e-12 + normals.Elem(i).Length());
  
	  for (i = 1; i <= nse; i++)
	  {
	  Element2d & el = SurfaceElement(i);
	  Vec3d nf = Cross (Vec3d (Point (el.PNum(1)), Point(el.PNum(2))),
	  Vec3d (Point (el.PNum(1)), Point(el.PNum(3))));
	  nf /= nf.Length();
	  Point3d c = Center (Point(el.PNum(1)),
	  Point(el.PNum(2)),
	  Point(el.PNum(3)));
  
	  for (j = 1; j <= 3; j++)
	  {
	  if (!linepoint.Test (el.PNum(j)))
	  {
	  double dist = Dist (c, Point(el.PNum(j)));
	  double dn = (nf - normals.Get(el.PNum(j))).Length();
  
	  RestrictLocalH (Point(el.PNum(j)), dist / (dn+1e-12) /elperr);
	  }
	  }
	  }
	  */
	}

	///
	public void CalcLocalHFromPointDistances(double grading)
	{
	  PrintMessage(3, "Calculating local h from point distances");

	  if (lochfunc == null)
	  {
		  Point3d pmin = new Point3d();
		  Point3d pmax = new Point3d();
		  GetBox(ref pmin, ref pmax);

		  // SetLocalH (pmin, pmax, mparam.grading);
	  SetLocalH(new netgen.Point3d(pmin), new netgen.Point3d(pmax), grading);
	  }

	  PointIndex i = new PointIndex();
	  PointIndex j = new PointIndex();
	  double hl;


	  for (i = PointIndex.BASE; i < GetNP() + PointIndex.BASE; i++)
	  {
		  for (j = i + 1; j < GetNP() + PointIndex.BASE; j++)
		  {
			  Point3d p1 = points[i];
			  Point3d p2 = points[j];
			  hl = netgen.GlobalMembers.Dist(p1,p2);
			  RestrictLocalH(p1, hl);
			  RestrictLocalH(p2, hl);
			  //cout << "restricted h at " << p1 << " and " << p2 << " to " << hl << endl;
		  }
	  }


	}

	///
	public void RestrictLocalH(resthtype rht, int nr, double loch)
	{
	  int i;
	  switch (rht)
	  {
		case resthtype.RESTRICTH_FACE:
		{
			for (i = 1; i <= GetNSE(); i++)
			{
				Element2d sel = SurfaceElement(i);
				if (sel.GetIndex() == nr)
				{
				  RestrictLocalH(resthtype.RESTRICTH_SURFACEELEMENT, i, loch);
				}
			}
			break;
		}
		case resthtype.RESTRICTH_EDGE:
		{
			for (i = 1; i <= GetNSeg(); i++)
			{
				Segment seg = LineSegment(i);
				if (seg.edgenr == nr)
				{
				  RestrictLocalH(resthtype.RESTRICTH_SEGMENT, i, loch);
				}
			}
			break;
		}
		case resthtype.RESTRICTH_POINT:
		{
			RestrictLocalH(new Point(nr), loch);
			break;
		}

		case resthtype.RESTRICTH_SURFACEELEMENT:
		{
			Element2d sel = SurfaceElement(nr);
			Point3d p = netgen.GlobalMembers.Center(new Point(sel.PNum(1)), new Point(sel.PNum(2)), new Point(sel.PNum(3)));
			RestrictLocalH(p, loch);
			break;
		}
		case resthtype.RESTRICTH_SEGMENT:
		{
			Segment seg = LineSegment(nr);
			RestrictLocalHLine(new Point(seg[0]), new Point(seg[1]), loch);
			break;
		}
	  }
	}

	///
	public void LoadLocalMeshSize(string meshsizefilename)
	{
	  // Philippose - 10/03/2009
	  // Improve error checking when loading and reading
	  // the local mesh size file

	  if (string.IsNullOrEmpty(meshsizefilename))
	  {
		  return;
	  }

	  ifstream msf = new ifstream(meshsizefilename);

	  // Philippose - 09/03/2009
	  // Adding print message information in case the specified
	  // does not exist, or does not load successfully due to
	  // other reasons such as access rights, etc...
	  if (msf == null)
	  {
		  PrintMessage(3, "Error loading mesh size file: ", meshsizefilename, "....","Skipping!");
		  return;
	  }

	  PrintMessage(3, "Load local mesh-size file: ", meshsizefilename);

	  int nmsp = 0;
	  int nmsl = 0;

	  msf >> nmsp;
	  if (!msf.good())
	  {
		throw new Exception("Mesh-size file error: No points found\n");
	  }

	  if (nmsp > 0)
	  {
		PrintMessage(4, "Number of mesh-size restriction points: ", nmsp);
	  }

	  for (int i = 0; i < nmsp; i++)
	  {
		  Point3d pi = new Point3d();
		  double hi;
		  msf >> pi.X() >> pi.Y() >> pi.Z();
		  msf >> hi;
		  if (!msf.good())
		  {
			throw new Exception("Mesh-size file error: Number of points don't match specified list size\n");
		  }
		  RestrictLocalH(pi, hi);
	  }

	  msf >> nmsl;
	  if (!msf.good())
	  {
		throw new Exception("Mesh-size file error: No line definitions found\n");
	  }

	  if (nmsl > 0)
	  {
		PrintMessage(4, "Number of mesh-size restriction lines: ", nmsl);
	  }

	  for (int i = 0; i < nmsl; i++)
	  {
		  Point3d p1 = new Point3d();
		  Point3d p2 = new Point3d();
		  double hi;
		  msf >> p1.X() >> p1.Y() >> p1.Z();
		  msf >> p2.X() >> p2.Y() >> p2.Z();
		  msf >> hi;
		  if (!msf.good())
		  {
			throw new Exception("Mesh-size file error: Number of line definitions don't match specified list size\n");
		  }
		  RestrictLocalHLine(p1, p2, hi);
	  }

	  msf.close();
	}

	///
	public void SetGlobalH(double h)
	{
	  hglob = h;
	}

	///
	public void SetMinimalH(double h)
	{
	  hmin = h;
	}

	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: double MaxHDomain(int dom) const
	public double MaxHDomain(int dom)
	{
	  if (maxhdomain.Size())
	  {
		return maxhdomain.Get(dom);
	  }
	  else
	  {
		return 1e10;
	  }
	}

	///
	public void SetMaxHDomain(Array<double> mhd)
	{
	  maxhdomain.SetSize(mhd.Size());
	  for (int i = 1; i <= mhd.Size(); i++)
	  {
		maxhdomain.Elem(i) = mhd.Get(i);
	  }
	}

	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: double GetH(const Point3d & p) const
	public double GetH(Point3d p)
	{
	  double hmin = hglob;
	  if (lochfunc != null)
	  {
		  double hl = lochfunc.GetH(new netgen.Point3d(p));
		  if (hl < hglob)
		  {
			hmin = hl;
		  }
	  }
	  return hmin;
	}

	///
	public double GetMinH(Point3d pmin, Point3d pmax)
	{
	  double hmin = hglob;
	  if (lochfunc != null)
	  {
		  double hl = lochfunc.GetMinH(new netgen.Point3d(pmin), new netgen.Point3d(pmax));
		  if (hl < hmin)
		  {
			hmin = hl;
		  }
	  }
	  return hmin;
	}

	///
	public bool HasLocalHFunction()
	{
		return lochfunc != null;
	}
	///
	public LocalH LocalHFunction()
	{
		return  lochfunc;
	}
	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool LocalHFunctionGenerated() const
	public bool LocalHFunctionGenerated()
	{
		return (lochfunc != null);
	}

	/// Find bounding box
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetBox(Point3d & pmin, Point3d & pmax, int dom = -1) const
	public void GetBox(ref Point3d pmin, ref Point3d pmax, int dom = -1)
	{
	  if (points.Size() == 0)
	  {
		  pmin = pmax = new Point3d(0, 0, 0);
		  return;
	  }

	  if (dom <= 0)
	  {
		  pmin = new Point3d(1e10, 1e10, 1e10);
		  pmax = new Point3d(-1e10, -1e10, -1e10);

		  for (PointIndex pi = points.Begin(); pi < points.End(); pi++)
		  {
			  pmin.SetToMin(this[pi]);
			  pmax.SetToMax(this[pi]);
		  }
	  }
	  else
	  {
		  int j;
		  int nse = GetNSE();
		  SurfaceElementIndex sei = new SurfaceElementIndex();

		  pmin = new Point3d(1e10, 1e10, 1e10);
		  pmax = new Point3d(-1e10, -1e10, -1e10);
		  for (sei = 0; sei < nse; sei++)
		  {
			  Element2d el = this[sei];
			  if (el.IsDeleted())
			  {
				  continue;
			  }

			  if (dom == -1 || el.GetIndex() == dom)
			  {
				  for (j = 0; j < 3; j++)
				  {
					  pmin.SetToMin(this[el[j]]);
					  pmax.SetToMax(this[el[j]]);
				  }
			  }
		  }
	  }

	  if (pmin.X() > 0.5e10)
	  {
		  pmin = pmax = new Point3d(0, 0, 0);
	  }
	}

	/// Find bounding box of points of typ ptyp or less
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetBox(Point3d & pmin, Point3d & pmax, POINTTYPE ptyp) const
	public void GetBox(ref Point3d pmin, ref Point3d pmax, POINTTYPE ptyp)
	{
	  if (points.Size() == 0)
	  {
		  pmin = pmax = new Point3d(0, 0, 0);
		  return;
	  }

	  pmin = new Point3d(1e10, 1e10, 1e10);
	  pmax = new Point3d(-1e10, -1e10, -1e10);

	  for (PointIndex pi = points.Begin(); pi < points.End(); pi++)
	  {
		if (points[pi].Type() <= ptyp)
		{
			pmin.SetToMin(this[pi]);
			pmax.SetToMax(this[pi]);
		}
	  }
	}

	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetNOpenElements() const
	public int GetNOpenElements()
	{
		return openelements.Size();
	}
	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const Element2d & OpenElement(int i) const
	public Element2d OpenElement(int i)
	{
		return openelements.Get(i);
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: auto & OpenElements() const
	public Array<Element2d> OpenElements()
	{
		return new Array<Element2d>(openelements);
	}

	/// are also quads open elements
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool HasOpenQuads() const
	public bool HasOpenQuads()
	{
	  int no = GetNOpenElements();
	  for (int i = 0; i < no; i++)
	  {
		if (openelements[i].GetNP() == 4)
		{
		  return true;
		}
	  }
	  return false;
	}

	/// split into connected pieces
	public void SplitIntoParts()
	{
	  int i;
	  int j;
	  int dom;
	  int ne = GetNE();
	  int np = GetNP();
	  int nse = GetNSE();

	  BitArray surfused = new BitArray(nse);
	  BitArray pused = new BitArray(np);

	  surfused.Clear();

	  dom = 0;

	  while (true)
	  {
		  int cntd = 1;

		  dom++;

		  pused.Clear();

		  int found = 0;
		  for (i = 1; i <= nse; i++)
		  {
			if (!surfused.Test(i))
			{
				SurfaceElement(i).SetIndex(dom);
				for (j = 1; j <= 3; j++)
				{
				  pused.Set(SurfaceElement(i).PNum(j));
				}
				found = 1;
				cntd = 1;
				surfused.Set(i);
				break;
			}
		  }

		  if (found == 0)
		  {
			break;
		  }

		  int change;
		  do
		  {
			  change = 0;
			  for (i = 1; i <= nse; i++)
			  {
				  int @is = 0;
				  int isnot = 0;
				  for (j = 1; j <= 3; j++)
				  {
					if (pused.Test(SurfaceElement(i).PNum(j)))
					{
					  @is = 1;
					}
					else
					{
					  isnot = 1;
					}
				  }

				  if (@is != 0 && isnot != 0)
				  {
					  change = 1;
					  for (j = 1; j <= 3; j++)
					  {
						pused.Set(SurfaceElement(i).PNum(j));
					  }
				  }

				  if (@is != 0)
				  {
					  if (!surfused.Test(i))
					  {
						  surfused.Set(i);
						  SurfaceElement(i).SetIndex(dom);
						  cntd++;
					  }
				  }
			  }


			  for (i = 1; i <= ne; i++)
			  {
				  int @is = 0;
				  int isnot = 0;
				  for (j = 1; j <= 4; j++)
				  {
					if (pused.Test(VolumeElement(i).PNum(j)))
					{
					  @is = 1;
					}
					else
					{
					  isnot = 1;
					}
				  }

				  if (@is != 0 && isnot != 0)
				  {
					  change = 1;
					  for (j = 1; j <= 4; j++)
					  {
						pused.Set(VolumeElement(i).PNum(j));
					  }
				  }

				  if (@is != 0)
				  {
					  VolumeElement(i).SetIndex(dom);
				  }
			  }
		  } while (change != 0);

		  PrintMessage(3, "domain ", dom, " has ", cntd, " surfaceelements");
	  }

	  /*
	    facedecoding.SetSize (dom);
	    for (i = 1; i <= dom; i++)
	    {
	    facedecoding.Elem(i).surfnr = 0;
	    facedecoding.Elem(i).domin = i;
	    facedecoding.Elem(i).domout = 0;
	    }
	  */
	  ClearFaceDescriptors();
	  for (i = 1; i <= dom; i++)
	  {
		AddFaceDescriptor(new FaceDescriptor(0, i, 0, 0));
	  }
	  CalcSurfacesOfNode();
	  timestamp = netgen.GlobalMembers.NextTimeStamp();
	}

	/// 
	public void SplitSeparatedFaces()
	{
	  PrintMessage(3, "SplitSeparateFaces");
	  int fdi;
	  int np = GetNP();

	  BitArray usedp = new BitArray(np);
	  Array<SurfaceElementIndex> els_of_face = new Array<SurfaceElementIndex>();

	  fdi = 1;
	  while (fdi <= GetNFD())
	  {
		  GetSurfaceElementsOfFace(fdi, els_of_face);

		  if (els_of_face.Size() == 0)
		  {
			  continue;
		  }

		  SurfaceElementIndex firstel = els_of_face[0];

		  usedp.Clear();
//C++ TO C# CONVERTER TODO TASK: The following line was determined to contain a copy constructor call - this should be verified and a copy constructor should be created:
//ORIGINAL LINE: for (int j = 1; j <= SurfaceElement(firstel).GetNP(); j++)
		  for (int j = 1; j <= SurfaceElement(new netgen.SurfaceElementIndex(firstel)).GetNP(); j++)
		  {
//C++ TO C# CONVERTER TODO TASK: The following line was determined to contain a copy constructor call - this should be verified and a copy constructor should be created:
//ORIGINAL LINE: usedp.Set(SurfaceElement(firstel).PNum(j));
			usedp.Set(SurfaceElement(new netgen.SurfaceElementIndex(firstel)).PNum(j));
		  }

		  bool changed;
		  do
		  {
			  changed = false;

			  for (int i = 0; i < els_of_face.Size(); i++)
			  {
				  Element2d el = SurfaceElement(els_of_face[i]);

				  bool has = false;
				  bool hasno = false;
				  for (int j = 0; j < el.GetNP(); j++)
				  {
					  if (usedp.Test(new netgen.Element2d(el[j])))
					  {
						has = true;
					  }
					  else
					  {
						hasno = true;
					  }
				  }

				  if (has && hasno)
				  {
					changed = true;
				  }

				  if (has)
				  {
					for (int j = 0; j < el.GetNP(); j++)
					{
					  usedp.Set(new netgen.Element2d(el[j]));
					}
				  }
			  }
		  } while (changed);

		  int nface = 0;
		  for (int i = 0; i < els_of_face.Size(); i++)
		  {
			  Element2d el = SurfaceElement(els_of_face[i]);

			  int hasno = 0;
			  for (int j = 1; j <= el.GetNP(); j++)
			  {
				if (!usedp.Test(el.PNum(j)))
				{
				  hasno = 1;
				}
			  }

			  if (hasno != 0)
			  {
				  if (nface == 0)
				  {
					  FaceDescriptor nfd = GetFaceDescriptor(fdi);
					  nface = AddFaceDescriptor(nfd);
				  }

				  el.SetIndex(nface);
			  }
		  }

		  // reconnect list
		  if (nface != 0)
		  {
			  facedecoding[nface-1].firstelement = -1;
			  facedecoding[fdi - 1].firstelement = -1;

			  for (int i = 0; i < els_of_face.Size(); i++)
			  {
				  int ind = SurfaceElement(els_of_face[i]).GetIndex();
				  SurfaceElement(els_of_face[i]).next = facedecoding[ind - 1].firstelement;
				  facedecoding[ind - 1].firstelement = els_of_face[i];
			  }
		  }

		  fdi++;
	  }


	  /*
	    fdi = 1;
	    while (fdi <= GetNFD())
	    {
	    int firstel = 0;
	    for (int i = 1; i <= GetNSE(); i++)
	    if (SurfaceElement(i).GetIndex() == fdi)
	    {
	    firstel = i;
	    break;
	    }
	    if (!firstel) continue;
  
	    usedp.Clear();
	    for (int j = 1; j <= SurfaceElement(firstel).GetNP(); j++)
	    usedp.Set (SurfaceElement(firstel).PNum(j));
  
	    int changed;
	    do
	    {
	    changed = 0;
	    for (int i = 1; i <= GetNSE(); i++)
	    {
	    const Element2d & el = SurfaceElement(i);
	    if (el.GetIndex() != fdi)
	    continue;
  
	    int has = 0;
	    int hasno = 0;
	    for (int j = 1; j <= el.GetNP(); j++)
	    {
	    if (usedp.Test(el.PNum(j)))
	    has = 1;
	    else
	    hasno = 1;
	    }
	    if (has && hasno)
	    changed = 1;
  
	    if (has)
	    for (int j = 1; j <= el.GetNP(); j++)
	    usedp.Set (el.PNum(j));
	    }
	    }
	    while (changed);
  
	    int nface = 0;
	    for (int i = 1; i <= GetNSE(); i++)
	    {
	    Element2d & el = SurfaceElement(i);
	    if (el.GetIndex() != fdi)
	    continue;
  
	    int hasno = 0;
	    for (int j = 1; j <= el.GetNP(); j++)
	    {
	    if (!usedp.Test(el.PNum(j)))
	    hasno = 1;
	    }
  
	    if (hasno)
	    {
	    if (!nface)
	    {
	    FaceDescriptor nfd = GetFaceDescriptor(fdi);
	    nface = AddFaceDescriptor (nfd);
	    }
  
	    el.SetIndex (nface);
	    }
	    }
	    fdi++;
	    }
	  */
	}

	/// Refines mesh and projects points to true surface
	// void Refine (int levels, const CSGeometry * geom);


//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool BoundaryEdge(PointIndex pi1, PointIndex pi2) const
	public bool BoundaryEdge(PointIndex pi1, PointIndex pi2)
	{
	  if (boundaryedges == null)
	  {
//C++ TO C# CONVERTER TODO TASK: There is no equivalent to 'const_cast' in C#:
	const_cast<Mesh>(this).BuildBoundaryEdges();
	  }

	  INDEX_2 i2 = new INDEX_2(pi1, pi2);
	  i2.Sort();
	  return boundaryedges.Used(i2);
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool IsSegment(PointIndex pi1, PointIndex pi2) const
	public bool IsSegment(PointIndex pi1, PointIndex pi2)
	{
	  INDEX_2 i2 = new INDEX_2(pi1, pi2);
	  i2.Sort();
	  return segmentht.Used(i2);
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: SegmentIndex SegmentNr(PointIndex pi1, PointIndex pi2) const
	public SegmentIndex SegmentNr(PointIndex pi1, PointIndex pi2)
	{
	  INDEX_2 i2 = new INDEX_2(pi1, pi2);
	  i2.Sort();
	  return segmentht.Get(i2);
	}


	/**
	   Remove unused points. etc.
	*/
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
//C++ TO C# CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
	Timer Compress_t("Mesh::Compress");

	public void Compress()
	{
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	  static Timer t("Mesh::Compress");
	  RegionTimer reg = new RegionTimer(Compress_t);

	  Array<PointIndex,PointIndex.BASE,PointIndex> op2np = new Array<PointIndex,PointIndex.BASE,PointIndex>(GetNP());
	  Array<MeshPoint> hpoints = new Array<MeshPoint>();
	  BitArrayChar<PointIndex.BASE> pused = new BitArrayChar<PointIndex.BASE>(GetNP());

	  /*
	    (*testout) << "volels: " << endl;
	    for (i = 1; i <= volelements.Size(); i++)
	    {
	    for (j = 1; j <= volelements.Get(i).GetNP(); j++)
	    (*testout) << volelements.Get(i).PNum(j) << " ";
	    (*testout) << endl;
	    }
	    (*testout) << "np: " << GetNP() << endl;
	  */

	  for (int i = 0; i < volelements.Size(); i++)
	  {
		if (volelements[i][0] <= PointIndex.BASE-1 || volelements[i].IsDeleted())
		{
			volelements.Delete(i);
			i--;
		}
	  }


	  for (int i = 0; i < surfelements.Size(); i++)
	  {
		if (surfelements[i].IsDeleted())
		{
			surfelements.Delete(i);
			i--;
		}
	  }

	  for (int i = 0; i < segments.Size(); i++)
	  {
		if (segments[i][0] <= PointIndex.BASE-1)
		{
			segments.Delete(i);
			i--;
		}
	  }

	  for (int i = 0; i < segments.Size(); i++)
	  {
		if (segments[i].edgenr < 0)
		{
			segments.Delete(i--);
		}
	  }

	  pused.Clear();
	  for (int i = 0; i < volelements.Size(); i++)
	  {
		  Element el = volelements[i];
		  for (int j = 0; j < el.GetNP(); j++)
		  {
			pused.Set(el[j]);
		  }
	  }

	  for (int i = 0; i < surfelements.Size(); i++)
	  {
		  Element2d el = surfelements[i];
		  for (int j = 0; j < el.GetNP(); j++)
		  {
			pused.Set(el[j]);
		  }
	  }

	  for (int i = 0; i < segments.Size(); i++)
	  {
		  Segment seg = segments[i];
		  pused.Set(seg[0]);
		  pused.Set(seg[1]);
	  }

	  for (int i = 0; i < openelements.Size(); i++)
	  {
		  Element2d el = openelements[i];
		  for (int j = 0; j < el.GetNP(); j++)
		  {
			pused.Set(el[j]);
		  }
	  }

	  for (int i = 0; i < lockedpoints.Size(); i++)
	  {
		pused.Set(lockedpoints[i]);
	  }


	  /*
	  // compress points doesn't work for identified points !
	  if (identifiedpoints)
	  {
	  for (i = 1; i <= identifiedpoints->GetNBags(); i++)
	  if (identifiedpoints->GetBagSize(i))
	  {
	  pused.Set ();
	  break;
	  }
	  }
	  */
	  //  pused.Set();


	  int npi = PointIndex.BASE-1;

	  for (PointIndex pi = points.Begin(); pi < points.End(); pi++)
	  {
		if (pused.Test(pi) != 0)
		{
			npi++;
			op2np[pi] = npi;
			hpoints.Append(points[pi]);
		}
		else
		{
		  op2np[pi] = -1;
		}
	  }



	  points.SetSize(0);
	  for (int i = 0; i < hpoints.Size(); i++)
	  {
		points.Append(hpoints[i]);
	  }

	  for (int i = 1; i <= volelements.Size(); i++)
	  {
		  Element el = VolumeElement(i);
		  for (int j = 0; j < el.GetNP(); j++)
		  {
			el[j] = op2np[el[j]];
		  }
	  }

	  for (int i = 1; i <= surfelements.Size(); i++)
	  {
		  Element2d el = SurfaceElement(i);
		  for (int j = 0; j < el.GetNP(); j++)
		  {
			el[j] = op2np[el[j]];
		  }
	  }

	  for (int i = 0; i < segments.Size(); i++)
	  {
		  Segment seg = segments[i];
		  seg[0] = op2np[seg[0]];
		  seg[1] = op2np[seg[1]];
	  }

	  for (int i = 1; i <= openelements.Size(); i++)
	  {
		  Element2d el = openelements.Elem(i);
		  for (int j = 0; j < el.GetNP(); j++)
		  {
			el[j] = op2np[el[j]];
		  }
	  }


	  for (int i = 0; i < lockedpoints.Size(); i++)
	  {
		lockedpoints[i] = op2np[lockedpoints[i]];
	  }
	  /*
	  for (int i = 0; i < facedecoding.Size(); i++)
	    facedecoding[i].firstelement = -1;
	  for (int i = surfelements.Size()-1; i >= 0; i--)
	    {
	      int ind = surfelements[i].GetIndex();
	      surfelements[i].next = facedecoding[ind-1].firstelement;
	      facedecoding[ind-1].firstelement = i;
	    }
	  */
	  RebuildSurfaceElementLists();
	  CalcSurfacesOfNode();


	  //  FindOpenElements();
	  timestamp = netgen.GlobalMembers.NextTimeStamp();
	}

	/// first vertex has lowest index
	public void OrderElements()
	{
	  foreach (var el in surfelements)
	  {
		  if (el.GetType() == ELEMENT_TYPE.TRIG)
		  {
			while (el[0] > el[1] || el[0] > el[2])
			{ // rotate element
				var hp = el[0];
				el[0] = el[1];
				el[1] = el[2];
				el[2] = hp;
				var hgi = el.GeomInfoPi(1);
				el.GeomInfoPi(1) = el.GeomInfoPi(2);
				el.GeomInfoPi(2) = el.GeomInfoPi(3);
				el.GeomInfoPi(3) = hgi;
			}
		  }
	  }

	  foreach (var el in volelements)
	  {
		if (el.GetType() == ELEMENT_TYPE.TET)
		{
			// lowest index first ...
			int mini = 0;
			for (int i = 1; i < 4; i++)
			{
			  if (el[i] < el[mini])
			  {
				  mini = i;
			  }
			}
			if (mini != 0)
			{ // swap 0 with mini, and the other two ...
				int i3 = -1;
				int i4 = -1;
				for (int i = 1; i < 4; i++)
				{
				  if (i != mini)
				  {
					  i4 = i3;
					  i3 = i;
				  }
				}
				swap(el[0], el[mini]);
				swap(el[i3], el[i4]);
			}

			while (el[1] > el[2] || el[1] > el[3])
			{ // rotate element to move second index to second position
				var hp = el[1];
				el[1] = el[2];
				el[2] = el[3];
				el[3] = hp;
			}
		}
	  }
	}

	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void Save(ostream & outfile) const
	public void Save(ostream outfile)
	{
	  int i;
	  int j;

	  double scale = 1; // globflags.GetNumFlag ("scale", 1);
	  int inverttets = 0; // globflags.GetDefineFlag ("inverttets");
	  int invertsurf = 0; // globflags.GetDefineFlag ("invertsurfacemesh");



	  outfile << "mesh3d" << "\n";

	  outfile << "dimension\n" << GetDimension() << "\n";

	  outfile << "geomtype\n" << (int)geomtype << "\n";


	  outfile << "\n";
	  outfile << "# surfnr    bcnr   domin  domout      np      p1      p2      p3" << "\n";


	  switch (geomtype)
	  {
		case GEOM_TYPE.GEOM_STL:
		  outfile << "surfaceelementsgi" << "\n";
		  break;
		case GEOM_TYPE.GEOM_OCC:
	case GEOM_TYPE.GEOM_ACIS:
		  outfile << "surfaceelementsuv" << "\n";
		  break;
		default:
		  outfile << "surfaceelements" << "\n";
		  break;
	  }

	  outfile << GetNSE() << "\n";

	  SurfaceElementIndex sei = new SurfaceElementIndex();
	  for (sei = 0; sei < GetNSE(); sei++)
	  {
		  if (this[sei].GetIndex())
		  {
			  outfile << " " << GetFaceDescriptor(this[sei].GetIndex()).SurfNr() + 1;
			  outfile << " " << GetFaceDescriptor(this[sei].GetIndex()).BCProperty();
			  outfile << " " << GetFaceDescriptor(this[sei].GetIndex()).DomainIn();
			  outfile << " " << GetFaceDescriptor(this[sei].GetIndex()).DomainOut();
		  }
		  else
		  {
			outfile << " 0 0 0";
		  }

		  Element2d sel = this[sei];
		  if (invertsurf != 0)
		  {
			sel.Invert();
		  }

		  outfile << " " << sel.GetNP();
		  for (j = 0; j < sel.GetNP(); j++)
		  {
			outfile << " " << sel[j];
		  }

		  switch (geomtype)
		  {
			case GEOM_TYPE.GEOM_STL:
			  for (j = 1; j <= sel.GetNP(); j++)
			  {
				outfile << " " << sel.GeomInfoPi(j).trignum;
			  }
			  break;
			case GEOM_TYPE.GEOM_OCC:
		case GEOM_TYPE.GEOM_ACIS:
			  for (j = 1; j <= sel.GetNP(); j++)
			  {
				  outfile << " " << sel.GeomInfoPi(j).u;
				  outfile << " " << sel.GeomInfoPi(j).v;
			  }
			  break;
			default:
			  ;
			  break;
		  }
		  outfile << "\n";
	  }

	  outfile << "\n" << "\n";
	  outfile << "#  matnr      np      p1      p2      p3      p4" << "\n";
	  outfile << "volumeelements" << "\n";
	  outfile << GetNE() << "\n";

	  for (ElementIndex ei = 0; ei < GetNE(); ei++)
	  {
		  outfile << this[ei].GetIndex();
		  outfile << " " << this[ei].GetNP();

		  Element el = this[ei];
		  if (inverttets != 0)
		  {
			  el.Invert();
		  }

		  for (j = 0; j < el.GetNP(); j++)
		  {
		outfile << " " << el[j];
		  }
		  outfile << "\n";
	  }


	  outfile << "\n" << "\n";
	  //     outfile << "   surf1   surf2      p1      p2" << "\n";
	  outfile << "# surfid  0   p1   p2   trignum1    trignum2   domin/surfnr1    domout/surfnr2   ednr1   dist1   ednr2   dist2 \n";
	  outfile << "edgesegmentsgi2" << "\n";
	  outfile << GetNSeg() << "\n";

	  for (i = 1; i <= GetNSeg(); i++)
	  {
		  Segment seg = LineSegment(i);
		  outfile.width(8);
		  outfile << seg.si; // 2D: bc number, 3D: wievielte Kante
		  outfile.width(8);
		  outfile << 0;
		  outfile.width(8);
		  outfile << seg[0];
		  outfile.width(8);
		  outfile << seg[1];
		  outfile << " ";
		  outfile.width(8);
		  outfile << seg.geominfo[0].trignum; // stl dreiecke
		  outfile << " ";
		  outfile.width(8);
		  outfile << seg.geominfo[1].trignum; // << endl;  // stl dreieck

		  if (dimension == 3)
		  {
			  outfile << " ";
			  outfile.width(8);
			  outfile << seg.surfnr1 + 1;
			  outfile << " ";
			  outfile.width(8);
			  outfile << seg.surfnr2 + 1;
		  }
		  else
		  {
			  outfile << " ";
			  outfile.width(8);
			  outfile << seg.domin;
			  outfile << " ";
			  outfile.width(8);
			  outfile << seg.domout;
		  }

		  outfile << " ";
		  outfile.width(8);
		  outfile << seg.edgenr;
		  outfile << " ";
		  outfile.width(12);
		  outfile.precision(16);
		  outfile << seg.epgeominfo[0].dist; // splineparameter (2D)
		  outfile << " ";
		  outfile.width(8);
		  outfile.precision(16);
		  outfile << seg.epgeominfo[1].edgenr; // geometry dependent
		  outfile << " ";
		  outfile.width(12);
		  outfile << seg.epgeominfo[1].dist;

		  outfile << "\n";
	  }


	  outfile << "\n" << "\n";
	  outfile << "#          X             Y             Z" << "\n";
	  outfile << "points" << "\n";
	  outfile << GetNP() << "\n";
	  outfile.precision(16);
	  outfile.setf(ios.@fixed, ios.floatfield);
	  outfile.setf(ios.showpoint);

	  PointIndex pi = new PointIndex();
	  for (pi = PointIndex.BASE; pi < GetNP() + PointIndex.BASE; pi++)
	  {
		  outfile.width(22);
		  outfile << this[pi](0) / scale << "  ";
		  outfile.width(22);
		  outfile << this[pi](1) / scale << "  ";
		  outfile.width(22);
		  outfile << this[pi](2) / scale << "\n";
	  }

	  if (ident.GetMaxNr() > 0)
	  {
		  outfile << "identifications\n";
		  Array<INDEX_2> identpairs = new Array<INDEX_2>();
		  int cnt = 0;
		  for (i = 1; i <= ident.GetMaxNr(); i++)
		  {
			  ident.GetPairs(i, identpairs);
			  cnt += identpairs.Size();
		  }
		  outfile << cnt << "\n";
		  for (i = 1; i <= ident.GetMaxNr(); i++)
		  {
			  ident.GetPairs(i, identpairs);
			  for (j = 1; j <= identpairs.Size(); j++)
			  {
				  outfile.width(8);
				  outfile << identpairs.Get(j).I1();
				  outfile.width(8);
				  outfile << identpairs.Get(j).I2();
				  outfile.width(8);
				  outfile << i << "\n";
			  }
		  }

		  outfile << "identificationtypes\n";
		  outfile << ident.GetMaxNr() << "\n";
		  for (i = 1; i <= ident.GetMaxNr(); i++)
		  {
			  int type = (int)ident.GetType(i);
			  outfile << " " << type;
		  }
		  outfile << "\n";
	  }

	  int cntmat = 0;
	  for (i = 1; i <= materials.Size(); i++)
	  {
		if (materials.Get(i) && (materials.Get(i).length()) != 0)
		{
		  cntmat++;
		}
	  }

	  if (cntmat != 0)
	  {
		  outfile << "materials" << "\n";
		  outfile << cntmat << "\n";
		  for (i = 1; i <= materials.Size(); i++)
		  {
			if (materials.Get(i) && (materials.Get(i).length()) != 0)
			{
			  outfile << i << " " << *materials.Get(i) << "\n";
			}
		  }
	  }


	  int cntbcnames = 0;
	  for (int ii = 0; ii < bcnames.Size(); ii++)
	  {
		if (bcnames[ii])
		{
			cntbcnames++;
		}
	  }

	  if (cntbcnames != 0)
	  {
		  outfile << "\n\nbcnames" << "\n" << bcnames.Size() << "\n";
		  for (i = 0; i < bcnames.Size(); i++)
		  {
			outfile << i + 1 << "\t" << GetBCName(i) << "\n";
		  }
		  outfile << "\n" << "\n";
	  }
	  int cntcd2names = 0;
	  for (int ii = 0; ii < cd2names.Size(); ii++)
	  {
		if (cd2names[ii])
		{
			cntcd2names++;
		}
	  }

	  if (cntcd2names != 0)
	  {
	  outfile << "\n\ncd2names" << "\n" << cd2names.Size() << "\n";
	  for (i = 0; i < cd2names.Size(); i++)
	  {
		outfile << i + 1 << "\t" << GetCD2Name(i) << "\n";
	  }
	  outfile << "\n" << "\n";
	  }

	  /*
	    if ( GetDimension() == 2 )
	    {
	    for (i = 1; i <= GetNSeg(); i++)
	    {
	    const Segment & seg = LineSegment (i);
	    if ( ! bcprops.Contains(seg.si) && seg.GetBCName() != "" )
	    {
	    bcprops.Append(seg.si);
	    cntbcnames++;
	    }
	    }
	    }
	    else
	    {
	    for (sei = 0; sei < GetNSE(); sei++)
	    {
	    if ((*this)[sei].GetIndex())
	    {
	    int bcp = GetFaceDescriptor((*this)[sei].GetIndex ()).BCProperty();
	    string name = GetFaceDescriptor((*this)[sei].GetIndex ()).BCName();
	    if ( !bcprops.Contains(bcp) &&
	    name != "" )
	    {
	    bcprops.Append(bcp);
	    cntbcnames++;
	    }
	    }
	    }
	    }
  
	    bcprops.SetSize(0);
	    if ( cntbcnames )
	    {
	    outfile << "\nbcnames" << endl << cntbcnames << endl;
	    if ( GetDimension() == 2 )
	    {
	    for (i = 1; i <= GetNSeg(); i++)
	    {
	    const Segment & seg = LineSegment (i);
	    if ( ! bcprops.Contains(seg.si) && seg.GetBCName() != "" )
	    {
	    bcprops.Append(seg.si);
	    outfile << seg.si << "\t" << seg.GetBCName() << endl;
	    }
	    }
	    }
	    else
	    {
	    for (sei = 0; sei < GetNSE(); sei++)
	    {
	    if ((*this)[sei].GetIndex())
	    {
	    int bcp = GetFaceDescriptor((*this)[sei].GetIndex ()).BCProperty();
	    string name = GetFaceDescriptor((*this)[sei].GetIndex ()).BCName();
	    if ( !bcprops.Contains(bcp) &&
	    name != "" )
	    {
	    bcprops.Append(bcp);
	    outfile << bcp << "\t" << name << endl;
	    }
	    }
	    }
	    }
	    outfile << endl << endl;
	    }
	  */

	  int cnt_sing = 0;
	  for (PointIndex pi = points.Begin(); pi < points.End(); pi++)
	  {
		if (this[pi].Singularity() >= 1.0)
		{
			cnt_sing++;
		}
	  }

	  if (cnt_sing != 0)
	  {
		  outfile << "singular_points" << "\n" << cnt_sing << "\n";
		  for (PointIndex pi = points.Begin(); pi < points.End(); pi++)
		  {
			if (this[pi].Singularity() >= 1.0)
			{
			  outfile << (int)pi << "\t" << this[pi].Singularity() << "\n";
			}
		  }
	  }

	  cnt_sing = 0;
	  for (SegmentIndex si = 0; si < GetNSeg(); si++)
	  {
		if (segments[si].singedge_left)
		{
			cnt_sing++;
		}
	  }
	  if (cnt_sing != 0)
	  {
		  outfile << "singular_edge_left" << "\n" << cnt_sing << "\n";
		  for (SegmentIndex si = 0; si < GetNSeg(); si++)
		  {
			if (segments[si].singedge_left)
			{
			  outfile << (int)si << "\t" << segments[si].singedge_left << "\n";
			}
		  }
	  }

	  cnt_sing = 0;
	  for (SegmentIndex si = 0; si < GetNSeg(); si++)
	  {
		if (segments[si].singedge_right)
		{
			cnt_sing++;
		}
	  }
	  if (cnt_sing != 0)
	  {
		  outfile << "singular_edge_right" << "\n" << cnt_sing << "\n";
		  for (SegmentIndex si = 0; si < GetNSeg(); si++)
		  {
			if (segments[si].singedge_right)
			{
			  outfile << (int)si << "\t" << segments[si].singedge_right << "\n";
			}
		  }
	  }


	  cnt_sing = 0;
	  for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
	  {
		if (GetFaceDescriptor(this[sei].GetIndex()).domin_singular)
		{
		  cnt_sing++;
		}
	  }

	  if (cnt_sing != 0)
	  {
		  outfile << "singular_face_inside" << "\n" << cnt_sing << "\n";
		  for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
		  {
			if (GetFaceDescriptor(this[sei].GetIndex()).domin_singular)
			{
			  outfile << (int)sei << "\t" << GetFaceDescriptor(this[sei].GetIndex()).domin_singular << "\n";
			}
		  }
	  }

	  cnt_sing = 0;
	  for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
	  {
		if (GetFaceDescriptor(this[sei].GetIndex()).domout_singular)
		{
			cnt_sing++;
		}
	  }
	  if (cnt_sing != 0)
	  {
		  outfile << "singular_face_outside" << "\n" << cnt_sing << "\n";
		  for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
		  {
			if (GetFaceDescriptor(this[sei].GetIndex()).domout_singular)
			{
			  outfile << (int)sei << "\t" << GetFaceDescriptor(this[sei].GetIndex()).domout_singular << "\n";
			}
		  }
	  }


	  // Philippose - 09/07/2009
	  // Add mesh face colours to Netgen Vol file format
	  // The colours are saved in RGB triplets
	  int cnt_facedesc = GetNFD();
	  if (cnt_facedesc != 0)
	  {
		 outfile << "\n" << "\n" << "#   Surfnr     Red     Green     Blue" << "\n";
		 outfile << "face_colours" << "\n" << cnt_facedesc << "\n";

		 outfile.precision(8);
		 outfile.setf(ios.@fixed, ios.floatfield);
		 outfile.setf(ios.showpoint);

		 for (i = 1; i <= cnt_facedesc; i++)
		 {
			outfile.width(8);
			outfile << GetFaceDescriptor(i).SurfNr() + 1 << " ";
			outfile.width(12);
			outfile << GetFaceDescriptor(i).SurfColour().X() << " ";
			outfile.width(12);
			outfile << GetFaceDescriptor(i).SurfColour().Y() << " ";
			outfile.width(12);
			outfile << GetFaceDescriptor(i).SurfColour().Z();
			outfile << "\n";
		 }
	  }

	  outfile << "\n" << "\n" << "endmesh" << "\n" << "\n";
	  if (geometry != null)
	  {
		geometry.SaveToMeshFile(outfile);
	  }
	}

	///
	public void Load(istream infile)
	{
	  if (!(infile.good()))
	  {
		  Console.Write("cannot load mesh");
		  Console.Write("\n");
		  throw new Exception("mesh file not found");
	  }

	  int rank = GetCommunicator().Rank();
	  int ntasks = GetCommunicator().Size();

	  string str = new string(new char[100]);
	  int i;
	  int n;

	  double scale = 1; // globflags.GetNumFlag ("scale", 1);
	  int inverttets = 0; // globflags.GetDefineFlag ("inverttets");
	  int invertsurf = 0; // globflags.GetDefineFlag ("invertsurfacemesh");


	  facedecoding.SetSize(0);

	  bool endmesh = false;


	  while (infile.good() && !endmesh)
	  {
		  infile >> str;
		  if (string.Compare(str, "dimension") == 0)
		  {
			  infile >> dimension;
		  }

		  if (string.Compare(str, "geomtype") == 0)
		  {
			  int hi;
			  infile >> hi;
			  geomtype = GEOM_TYPE(hi);
		  }


		  if (string.Compare(str, "surfaceelements") == 0 || string.Compare(str, "surfaceelementsgi") == 0 || string.Compare(str, "surfaceelementsuv") == 0)
		  {
			  infile >> n;
			  PrintMessage(3, n, " surface elements");

		  bool geominfo = string.Compare(str, "surfaceelementsgi") == 0;
		  bool uv = string.Compare(str, "surfaceelementsuv") == 0;


			  for (i = 1; i <= n; i++)
			  {
				  int surfnr;
				  int bcp;
				  int domin;
				  int domout;
				  int nep;
				  int faceind = 0;

				  infile >> surfnr >> bcp >> domin >> domout;
				  surfnr--;

		  bool invert_el = false;
		  /*
		  if (domin == 0)
		    {
		      invert_el = true;
		      Swap (domin, domout);
		    }
		  */

				  for (int j = 1; j <= facedecoding.Size(); j++)
				  {
					if (GetFaceDescriptor(j).SurfNr() == surfnr && GetFaceDescriptor(j).BCProperty() == bcp && GetFaceDescriptor(j).DomainIn() == domin && GetFaceDescriptor(j).DomainOut() == domout)
					{
					  faceind = j;
					}
				  }

		  // if (facedecoding.Size()) faceind = 1;   // for timing

				  if (faceind == 0)
				  {
					  faceind = AddFaceDescriptor(new FaceDescriptor(surfnr, domin, domout, 0));
					  GetFaceDescriptor(faceind).SetBCProperty(bcp);
				  }

				  infile >> nep;
				  if (nep == 0)
				  {
					  nep = 3;
				  }

				  Element2d tri = new Element2d(nep);
				  tri.SetIndex(faceind);

				  for (int j = 1; j <= nep; j++)
				  {
					infile >> tri.PNum(j);
				  }

				  if (geominfo)
				  {
					for (int j = 1; j <= nep; j++)
					{
					  infile >> tri.GeomInfoPi(j).trignum;
					}
				  }

				  if (uv)
				  {
					for (int j = 1; j <= nep; j++)
					{
					  infile >> tri.GeomInfoPi(j).u >> tri.GeomInfoPi(j).v;
					}
				  }

				  if (invertsurf != 0)
				  {
					  tri.Invert();
				  }
		  if (invert_el)
		  {
			  tri.Invert();
		  }

		  AddSurfaceElement(tri);
			  }
		  }

		  if (string.Compare(str, "volumeelements") == 0)
		  {
			  infile >> n;
			  PrintMessage(3, n, " volume elements");
			  for (i = 1; i <= n; i++)
			  {
				  Element el = new Element(ELEMENT_TYPE.TET);
				  int hi;
				  int nep;
				  infile >> hi;
				  if (hi == 0)
				  {
					  hi = 1;
				  }
				  el.SetIndex(hi);
				  infile >> nep;
				  el.SetNP(nep);
				  el.SetCurved(nep != 4);
				  for (int j = 0; j < nep; j++)
				  {
					infile >> (int)(el[j]);
				  }

				  if (inverttets != 0)
				  {
					el.Invert();
				  }

		  AddVolumeElement(el);
			  }
		  }


		  if (string.Compare(str, "edgesegments") == 0)
		  {
			  infile >> n;
			  for (i = 1; i <= n; i++)
			  {
				  Segment seg = new Segment();
				  int hi;
				  infile >> seg.si >> hi >> seg[0] >> seg[1];
				  AddSegment(seg);
			  }
		  }



		  if (string.Compare(str, "edgesegmentsgi") == 0)
		  {
			  infile >> n;
			  for (i = 1; i <= n; i++)
			  {
				  Segment seg = new Segment();
				  int hi;
				  infile >> seg.si >> hi >> seg[0] >> seg[1] >> seg.geominfo[0].trignum >> seg.geominfo[1].trignum;
				  AddSegment(seg);
			  }
		  }

		  if (string.Compare(str, "edgesegmentsgi2") == 0)
		  {
			  int a;
			  infile >> a;
			  n = a;

			  PrintMessage(3, n, " curve elements");

			  for (i = 1; i <= n; i++)
			  {
				  Segment seg = new Segment();
				  int hi;
				  infile >> seg.si >> hi >> seg[0] >> seg[1] >> seg.geominfo[0].trignum >> seg.geominfo[1].trignum >> seg.surfnr1 >> seg.surfnr2 >> seg.edgenr >> seg.epgeominfo[0].dist >> seg.epgeominfo[1].edgenr >> seg.epgeominfo[1].dist;

				  seg.epgeominfo[0].edgenr = seg.epgeominfo[1].edgenr;

				  seg.domin = seg.surfnr1;
				  seg.domout = seg.surfnr2;

				  seg.surfnr1--;
				  seg.surfnr2--;

				  AddSegment(seg);
			  }
		  }

		  if (string.Compare(str, "points") == 0)
		  {
			  infile >> n;
			  PrintMessage(3, n, " points");
			  for (i = 1; i <= n; i++)
			  {
				  Point3d p = new Point3d();
				  infile >> p.X() >> p.Y() >> p.Z();
				  p.X() *= scale;
				  p.Y() *= scale;
				  p.Z() *= scale;
				  AddPoint(p);
			  }
		  PrintMessage(3, n, " points done");
		  }

		  if (string.Compare(str, "identifications") == 0)
		  {
			  infile >> n;
			  PrintMessage(3, n, " identifications");
			  for (i = 1; i <= n; i++)
			  {
				  PointIndex pi1 = new PointIndex();
				  PointIndex pi2 = new PointIndex();
				  int ind;
				  infile >> pi1 >> pi2 >> ind;
				  ident.Add(new netgen.PointIndex(pi1), new netgen.PointIndex(pi2), ind);
			  }
		  }

		  if (string.Compare(str, "identificationtypes") == 0)
		  {
			  infile >> n;
			  PrintMessage(3, n, " identificationtypes");
			  for (i = 1; i <= n; i++)
			  {
				  int type;
				  infile >> type;
				  ident.SetType(i, Identifications.ID_TYPE(type));
			  }
		  }

		  if (string.Compare(str, "materials") == 0)
		  {
			  infile >> n;
			  for (i = 1; i <= n; i++)
			  {
				  int nr;
				  string mat;
				  infile >> nr >> mat;
				  SetMaterial(nr, mat);
			  }
		  }

		  if (string.Compare(str, "bcnames") == 0)
		  {
			  infile >> n;
			  Array<int,0> bcnrs = new Array<int,0>(n);
			  SetNBCNames(n);
			  for (i = 1; i <= n; i++)
			  {
				  string nextbcname;
				  infile >> bcnrs[i - 1] >> nextbcname;
				  bcnames[bcnrs[i - 1] - 1] = new string(nextbcname);
			  }

			  if (GetDimension() == 2)
			  {
				  for (i = 1; i <= GetNSeg(); i++)
				  {
					  Segment seg = LineSegment(i);
					  if (seg.si <= n)
					  {
						seg.SetBCName(bcnames[seg.si - 1]);
					  }
					  else
					  {
						seg.SetBCName(null);
					  }
				  }
			  }
			  else
			  {
				  for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
				  {
					  if (this[sei].GetIndex())
					  {
						  int bcp = GetFaceDescriptor(this[sei].GetIndex()).BCProperty();
						  if (bcp <= n)
						  {
							GetFaceDescriptor(this[sei].GetIndex()).SetBCName(bcnames[bcp - 1]);
						  }
						  else
						  {
							GetFaceDescriptor(this[sei].GetIndex()).SetBCName(null);
						  }

					  }
				  }

			  }
		  }

	  if (string.Compare(str, "cd2names") == 0)
	  {
		  infile >> n;
		  Array<int,0> cd2nrs = new Array<int,0>(n);
		  SetNCD2Names(n);
		  for (i = 1; i <= n; i++)
		  {
		  string nextcd2name;
		  infile >> cd2nrs[i - 1] >> nextcd2name;
		  cd2names[cd2nrs[i - 1] - 1] = new string(nextcd2name);
		  }
		  if (GetDimension() == 2)
		  {
		  throw new Exception("co dim 2 elements not implemented for dimension 2");
		  }
		  else
		  {
		  for (i = 1; i <= GetNSeg(); i++)
		  {
			  Segment seg = LineSegment(i);
			  if (seg.edgenr <= n)
			  {
				seg.SetBCName(cd2names[seg.edgenr - 1]);
			  }
			  else
			  {
				seg.SetBCName(null);
			  }
		  }
		  }
	  }

		  if (string.Compare(str, "singular_points") == 0)
		  {
			  infile >> n;
			  for (i = 1; i <= n; i++)
			  {
				  PointIndex pi = new PointIndex();
				  double s;
				  infile >> pi;
				  infile >> s;
				  this[pi].Singularity(s);
			  }
		  }

		  if (string.Compare(str, "singular_edge_left") == 0)
		  {
			  infile >> n;
			  for (i = 1; i <= n; i++)
			  {
				  SegmentIndex si = new SegmentIndex();
				  double s;
				  infile >> si;
				  infile >> s;
				  this[si].singedge_left = s;
			  }
		  }
		  if (string.Compare(str, "singular_edge_right") == 0)
		  {
			  infile >> n;
			  for (i = 1; i <= n; i++)
			  {
				  SegmentIndex si = new SegmentIndex();
				  double s;
				  infile >> si;
				  infile >> s;
				  this[si].singedge_right = s;
			  }
		  }

		  if (string.Compare(str, "singular_face_inside") == 0)
		  {
			  infile >> n;
			  for (i = 1; i <= n; i++)
			  {
				  SurfaceElementIndex sei = new SurfaceElementIndex();
				  double s;
				  infile >> sei;
				  infile >> s;
				  GetFaceDescriptor(this[sei].GetIndex()).domin_singular = s;
			  }
		  }

		  if (string.Compare(str, "singular_face_outside") == 0)
		  {
			  infile >> n;
			  for (i = 1; i <= n; i++)
			  {
				  SurfaceElementIndex sei = new SurfaceElementIndex();
				  double s;
				  infile >> sei;
				  infile >> s;
				  GetFaceDescriptor(this[sei].GetIndex()).domout_singular = s;
			  }
		  }

		  // Philippose - 09/07/2009
		  // Add mesh face colours to Netgen Vol file format
		  // The colours are read in as RGB triplets
		  if (string.Compare(str, "face_colours") == 0)
		  {
			 int cnt_facedesc = GetNFD();
			 infile >> n;
			 if (n == cnt_facedesc)
			 {
				for (i = 1; i <= n; i++)
				{
				   int surfnr = 0;
				   Vec3d surfcolour = new Vec3d(0.0, 1.0, 0.0);

				   infile >> surfnr >> surfcolour.X() >> surfcolour.Y() >> surfcolour.Z();

				   surfnr--;

				   if (surfnr > 0)
				   {
					  for (int facedesc = 1; facedesc <= cnt_facedesc; facedesc++)
					  {
						 if (surfnr == GetFaceDescriptor(facedesc).SurfNr())
						 {
							GetFaceDescriptor(facedesc).SetSurfColour(new netgen.Vec3d(surfcolour));
						 }
					  }
				   }
				}
			 }
		  }

		  if (string.Compare(str, "endmesh") == 0)
		  {
			endmesh = true;
		  }



		  str = "";
	  }




	  CalcSurfacesOfNode();

	  if (ntasks == 1) // sequential run only
	  {
	  topology.Update();
	  clusters.Update();
	  }

	  SetNextMajorTimeStamp();
	  //  PrintMemInfo (cout);
	}

	///
	public void Merge(istream infile, int surfindex_offset = 0)
	{
	  string str = new string(new char[100]);
	  int i;
	  int n;


	  int inverttets = 0; // globflags.GetDefineFlag ("inverttets");

	  int oldnp = GetNP();
	  int oldne = GetNSeg();
	  int oldnd = GetNDomains();

	  for (SurfaceElementIndex si = 0; si < GetNSE(); si++)
	  {
		for (int j = 1; j <= this[si].GetNP(); j++)
		{
			this[si].GeomInfoPi(j).trignum = -1;
		}
	  }

	  int max_surfnr = 0;
	  for (i = 1; i <= GetNFD(); i++)
	  {
		max_surfnr = netgen.GlobalMembers.max2(max_surfnr, GetFaceDescriptor(i).SurfNr());
	  }
	  max_surfnr++;

	  if (max_surfnr < surfindex_offset)
	  {
		  max_surfnr = surfindex_offset;
	  }


	  bool endmesh = false;

	  while (infile.good() && !endmesh)
	  {
		  infile >> str;

		  if (string.Compare(str, "surfaceelementsgi") == 0 || string.Compare(str, "surfaceelements") == 0)
		  {
			  infile >> n;
			  PrintMessage(3, n, " surface elements");
			  for (i = 1; i <= n; i++)
			  {
				  int j;
				  int surfnr;
				  int bcp;
				  int domin;
				  int domout;
				  int nep;
				  int faceind = 0;
				  infile >> surfnr >> bcp >> domin >> domout;

				  surfnr--;

				  if (domin > 0)
				  {
					  domin += oldnd;
				  }
				  if (domout > 0)
				  {
					  domout += oldnd;
				  }
				  surfnr += max_surfnr;


				  for (j = 1; j <= facedecoding.Size(); j++)
				  {
					if (GetFaceDescriptor(j).SurfNr() == surfnr && GetFaceDescriptor(j).BCProperty() == bcp && GetFaceDescriptor(j).DomainIn() == domin && GetFaceDescriptor(j).DomainOut() == domout)
					{
					  faceind = j;
					}
				  }

				  if (faceind == 0)
				  {
					  faceind = AddFaceDescriptor(new FaceDescriptor(surfnr, domin, domout, 0));
					  if (GetDimension() == 2)
					  {
						  bcp++;
					  }
					  GetFaceDescriptor(faceind).SetBCProperty(bcp);
				  }

				  infile >> nep;
				  if (nep == 0)
				  {
					  nep = 3;
				  }

				  Element2d tri = new Element2d(nep);
				  tri.SetIndex(faceind);

				  for (j = 1; j <= nep; j++)
				  {
					  infile >> tri.PNum(j);
					  tri.PNum(j) = tri.PNum(j) + oldnp;
				  }


				  if (string.Compare(str, "surfaceelementsgi") == 0)
				  {
					for (j = 1; j <= nep; j++)
					{
						infile >> tri.GeomInfoPi(j).trignum;
						tri.GeomInfoPi(j).trignum = -1;
					}
				  }

				  AddSurfaceElement(tri);
			  }
		  }


		  if (string.Compare(str, "edgesegments") == 0)
		  {
			  infile >> n;
			  for (i = 1; i <= n; i++)
			  {
				  Segment seg = new Segment();
				  int hi;
				  infile >> seg.si >> hi >> seg[0] >> seg[1];
				  seg[0] = seg[0] + oldnp;
				  seg[1] = seg[1] + oldnp;
				  AddSegment(seg);
			  }
		  }



		  if (string.Compare(str, "edgesegmentsgi") == 0)
		  {
			  infile >> n;
			  for (i = 1; i <= n; i++)
			  {
				  Segment seg = new Segment();
				  int hi;
				  infile >> seg.si >> hi >> seg[0] >> seg[1] >> seg.geominfo[0].trignum >> seg.geominfo[1].trignum;
				  seg[0] = seg[0] + oldnp;
				  seg[1] = seg[1] + oldnp;
				  AddSegment(seg);
			  }
		  }
		  if (string.Compare(str, "edgesegmentsgi2") == 0)
		  {
			  infile >> n;
			  PrintMessage(3, n, " curve elements");

			  for (i = 1; i <= n; i++)
			  {
				  Segment seg = new Segment();
				  int hi;
				  infile >> seg.si >> hi >> seg[0] >> seg[1] >> seg.geominfo[0].trignum >> seg.geominfo[1].trignum >> seg.surfnr1 >> seg.surfnr2 >> seg.edgenr >> seg.epgeominfo[0].dist >> seg.epgeominfo[1].edgenr >> seg.epgeominfo[1].dist;
				  seg.epgeominfo[0].edgenr = seg.epgeominfo[1].edgenr;

				  seg.surfnr1--;
				  seg.surfnr2--;

				  if (seg.surfnr1 >= 0)
				  {
					  seg.surfnr1 = seg.surfnr1 + max_surfnr;
				  }
				  if (seg.surfnr2 >= 0)
				  {
					  seg.surfnr2 = seg.surfnr2 + max_surfnr;
				  }
				  seg[0] = seg[0] + oldnp;
				  seg[1] = seg[1] + oldnp;
		  *testout << "old edgenr: " << seg.edgenr << "\n";
				  seg.edgenr = seg.edgenr + oldne;
		  *testout << "new edgenr: " << seg.edgenr << "\n";
				  seg.epgeominfo[1].edgenr = seg.epgeominfo[1].edgenr + oldne;

				  AddSegment(seg);
			  }
		  }

		  if (string.Compare(str, "volumeelements") == 0)
		  {
			  infile >> n;
			  PrintMessage(3, n, " volume elements");
			  for (i = 1; i <= n; i++)
			  {
				  Element el = new Element(ELEMENT_TYPE.TET);
				  int hi;
				  int nep;
				  infile >> hi;
				  if (hi == 0)
				  {
					  hi = 1;
				  }
				  el.SetIndex(hi + oldnd);
				  infile >> nep;
				  el.SetNP(nep);

				  for (int j = 0; j < nep; j++)
				  {
					  infile >> (int)(el[j]);
					  el[j] = el[j] + oldnp;
				  }

				  if (inverttets != 0)
				  {
					el.Invert();
				  }

				  AddVolumeElement(el);
			  }
		  }


		  if (string.Compare(str, "points") == 0)
		  {
			  infile >> n;
			  PrintMessage(3, n, " points");
			  for (i = 1; i <= n; i++)
			  {
				  Point3d p = new Point3d();
				  infile >> p.X() >> p.Y() >> p.Z();
				  AddPoint(p);
			  }
		  }


		  if (string.Compare(str, "endmesh") == 0)
		  {
			  endmesh = true;
		  }


		  if (string.Compare(str, "materials") == 0)
		  {
			  infile >> n;
			  for (i = 1; i <= n; i++)
			  {
				  int nr;
				  string mat;
				  infile >> nr >> mat;
				  SetMaterial(nr + oldnd, mat);
			  }
		  }


		  str = "";
	  }

	  CalcSurfacesOfNode();

	  topology.Update();
	  clusters.Update();

	  SetNextMajorTimeStamp();
	}

	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void Save(const string & filename) const
	public void Save(string filename)
	{
	  ostream outfile;
	  if (filename.IndexOf(".vol.gz") != -1)
	  {
		outfile = new ogzstream(filename);
	  }
	  else if (filename.IndexOf(".vol") != -1)
	  {
		outfile = new ofstream(filename);
	  }
	  else
	  {
		outfile = new ogzstream((filename + ".vol.gz").c_str());
	  }

	  Save(outfile);
	  outfile = null;
	}

	///
	public void Load(string filename)
	{
	  Console.Write("filename = ");
	  Console.Write(filename);
	  Console.Write("\n");
	  istream infile = null;

	  if (filename.IndexOf(".vol.gz") != -1)
	  {
		infile = new igzstream(filename);
	  }
	  else
	  {
		infile = new ifstream(filename);
	  }

	  // ifstream infile(filename.c_str());
	  if (!(infile.good()))
	  {
		throw new Exception("mesh file not found");
	  }

	  Load(infile);
	  infile = null;
	}

	///
	public void Merge(string filename, int surfindex_offset = 0)
	{
	  ifstream infile = new ifstream(filename);
	  if (!infile.good())
	  {
		throw new Exception("mesh file not found");
	  }

	  Merge(infile,surfindex_offset);

	}


	public void DoArchive(Archive archive)
	{
	  archive dimension;
	  archive points;
	  archive surfelements;
	  archive volelements;
	  archive segments;
	  archive facedecoding;
	  archive & materials & bcnames & cd2names & cd3names;

	  archive & *ident;

	  archive.Shallow(geometry);
	  archive & *curvedelems;

	  if (archive.Input())
	  {
	  int rank = GetCommunicator().Rank();
	  int ntasks = GetCommunicator().Size();

		  RebuildSurfaceElementLists();

		  CalcSurfacesOfNode();
		  if (ntasks == 1) // sequential run only
		  {
			  topology.Update();
			  clusters.Update();
		  }
		  SetNextMajorTimeStamp();
	  }
	}

	///
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//	DLL_HEADER void ImproveMesh(MeshingParameters mp, OPTIMIZEGOAL goal = OPTIMIZEGOAL.OPT_QUALITY);

	///
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//	void ImproveMeshJacobian(MeshingParameters mp, OPTIMIZEGOAL goal = OPTIMIZEGOAL.OPT_QUALITY, BitArray usepoint = null);
	///
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//	void ImproveMeshJacobianOnSurface(MeshingParameters mp, BitArray usepoint, Array< Vec<3> > nv, OPTIMIZEGOAL goal = OPTIMIZEGOAL.OPT_QUALITY, Array< Array<int,PointIndex::BASE> > idmaps = null);
	/**
	   free nodes in environment of openelements 
	   for optimiztion
	*/
	public void FreeOpenElementsEnvironment(int layers)
	{
	  int i;
	  int j;
	  int k;
	  PointIndex pi = new PointIndex();
	  const int large = 9999;
	  Array<int,PointIndex.BASE> dist = new Array<int,PointIndex.BASE>(GetNP());

	  dist = large;

	  for (int i = 1; i <= GetNOpenElements(); i++)
	  {
		  Element2d face = OpenElement(i);
		  for (j = 0; j < face.GetNP(); j++)
		  {
			dist[face[j]] = 1;
		  }
	  }

	  for (k = 1; k <= layers; k++)
	  {
		for (i = 1; i <= GetNE(); i++)
		{
			Element el = VolumeElement(i);
			if (el[0] == -1 || el.IsDeleted())
			{
				continue;
			}

			int elmin = large;
			for (j = 0; j < el.GetNP(); j++)
			{
			  if (dist[el[j]] < elmin)
			  {
				elmin = dist[el[j]];
			  }
			}

			if (elmin < large)
			{
				for (j = 0; j < el.GetNP(); j++)
				{
				  if (dist[el[j]] > elmin + 1)
				  {
					dist[el[j]] = elmin + 1;
				  }
				}
			}
		}
	  }

	  int cntfree = 0;
	  for (i = 1; i <= GetNE(); i++)
	  {
		  Element el = VolumeElement(i);
		  if (el[0] == -1 || el.IsDeleted())
		  {
			  continue;
		  }

		  int elmin = large;
		  for (j = 0; j < el.GetNP(); j++)
		  {
			if (dist[el[j]] < elmin)
			{
			  elmin = dist[el[j]];
			}
		  }

		  el.flags.@fixed = elmin > layers;
		  // eltyps.Elem(i) = (elmin <= layers) ?
		  // FREEELEMENT : FIXEDELEMENT;
		  if (elmin <= layers)
		  {
			cntfree++;
		  }
	  }

	  PrintMessage(5, "free: ", cntfree, ", fixed: ", GetNE() - cntfree);
	  (*testout) << "free: " << cntfree << ", fixed: " << GetNE() - cntfree << "\n";

	  for (pi = PointIndex.BASE; pi < GetNP() + PointIndex.BASE; pi++)
	  {
		  if (dist[pi] > layers + 1)
		  {
			points[pi].SetType(POINTTYPE.FIXEDPOINT);
		  }
	  }
	}

	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool LegalTet(Element & el) const
	public bool LegalTet(Element el)
	{
	  if (el.IllegalValid() != 0)
	  {
	return el.Illegal() == 0;
	  }
	  return LegalTet2(el);
	}
	///

	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool LegalTet2(Element & el) const
	public bool LegalTet2(Element el)
	{
	  // static int timer1 = NgProfiler::CreateTimer ("Legaltet2");

	  // Test, whether 4 points have a common surface plus
	  // at least 4 edges at the boundary

	  if (boundaryedges == null)
	  {
//C++ TO C# CONVERTER TODO TASK: There is no equivalent to 'const_cast' in C#:
		const_cast<Mesh>(this).BuildBoundaryEdges();
	  }


	  // non-tets are always legal
	  if (el.GetType() != ELEMENT_TYPE.TET)
	  {
		  el.SetLegal(1);
		  return true;
	  }

	  POINTTYPE[] pointtype = new POINTTYPE[4];
	  for (int i = 0; i < 4; i++)
	  {
		pointtype[i] = this[el[i]].Type();
	  }



	  // element has at least 2 inner points ---> legal
	  int cnti = 0;
	  for (int j = 0; j < 4; j++)
	  {
		if (pointtype[j] == POINTTYPE.INNERPOINT)
		{
			cnti++;
			if (cnti >= 2)
			{
				el.SetLegal(1);
				return true;
			}
		}
	  }



	  // which faces are boundary faces ?
	  int[] bface = new int[4];
	  for (int i = 0; i < 4; i++)
	  {
		  bface[i] = surfelementht.Used(INDEX_3.Sort(new netgen.Element(el[GlobalMembers.gftetfacesa[i][0]]), new netgen.Element(el[GlobalMembers.gftetfacesa[i][1]]), new netgen.Element(el[GlobalMembers.gftetfacesa[i][2]])));
	  }

	  int[][] bedge = RectangularArrays.RectangularIntArray(4, 4);
	  int[][] segedge = RectangularArrays.RectangularIntArray(4, 4);
	  int[][] pi3map =
	  {
		  new int[] {-1, 2, 1, 1},
		  new int[] {2, -1, 0, 0},
		  new int[] {1, 0, -1, 0},
		  new int[] {1, 0, 0, -1}
	  };

	  int[][] pi4map =
	  {
		  new int[] {-1, 3, 3, 2},
		  new int[] {3, -1, 3, 2},
		  new int[] {3, 3, -1, 1},
		  new int[] {2, 2, 1, -1}
	  };


	  for (int i = 0; i < 4; i++)
	  {
		for (int j = 0; j < i; j++)
		{
			bool sege = false;
			bool be = false;

			int pos = boundaryedges.Position0(INDEX_2.Sort(new netgen.Element(el[i]), new netgen.Element(el[j])));
			if (pos != -1)
			{
				be = true;
				if (boundaryedges.GetData0(pos) == 2)
				{
				  sege = true;
				}
			}

			segedge[j][i] = segedge[i][j] = sege;
			bedge[j][i] = bedge[i][j] = be;
		}
	  }

	  // two boundary faces and no edge is illegal
	  for (int i = 0; i < 3; i++)
	  {
		for (int j = i + 1; j < 4; j++)
		{
			if (bface[i] != 0 && bface[j] != 0)
			{
			  if (segedge[pi3map[i][j]][pi4map[i][j]] == 0)
			  {
				  // 2 boundary faces withoud edge in between
				  el.SetLegal(0);
				  return false;
			  }
			}
		}
	  }

	  // three boundary edges meeting in a Surface point
	  for (int i = 0; i < 4; i++)
	  {
		  if (pointtype[i] == POINTTYPE.SURFACEPOINT)
		  {
			  bool alledges = true;
			  for (int j = 0; j < 4; j++)
			  {
				if (j != i && bedge[i][j] == 0)
				{
					alledges = false;
					break;
				}
			  }
			  if (alledges)
			  {
				  // cout << "tet illegal due to unmarked node" << endl;
				  el.SetLegal(0);
				  return false;
			  }
		  }
	  }



	  for (int fnr = 0; fnr < 4; fnr++)
	  {
		if (bface[fnr] == 0)
		{
		  for (int i = 0; i < 4; i++)
		  {
			if (i != fnr)
			{
				int pi1 = pi3map[i][fnr];
				int pi2 = pi4map[i][fnr];

				if (pointtype[i] == POINTTYPE.SURFACEPOINT)
				{
					// two connected edges on surface, but no face
					if (bedge[i][pi1] != 0 && bedge[i][pi2] != 0)
					{
						el.SetLegal(0);
						return false;
					}
				}

				if (pointtype[i] == POINTTYPE.EDGEPOINT)
				{
					// connected surface edge and edge edge, but no face
					if ((bedge[i][pi1] && segedge[i][pi2]) || (bedge[i][pi2] && segedge[i][pi1]))
					{
						el.SetLegal(0);
						return false;
					}
				}

			}
		  }
		}
	  }


	  el.SetLegal(1);
	  return true;

	}


	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool LegalTrig(const Element2d & el) const
	public bool LegalTrig(Element2d el)
	{
	  return true;
	  if (true) // needed for old, simple hp-refinement
	  {
		  // trigs with 2 or more segments are illegal
		  int i;
		  int nseg = 0;

		  if (segmentht == null)
		  {
			  cerr << "no segmentht allocated" << "\n";
			  return false;
		  }

		  //      Point3d cp(0.5, 0.5, 0.5);
		  for (i = 1; i <= 3; i++)
		  {
			  INDEX_2 i2 = new INDEX_2(el.PNumMod(i), el.PNumMod(i + 1));
			  i2.Sort();
			  if (segmentht.Used(i2))
			  {
				nseg++;
			  }
		  }
		  if (nseg >= 2)
		  {
			return false;
		  }
	  }
	  return true;
	}

	/**
	   if values non-null, return values in 4-double array:
	   triangle angles min/max, tetangles min/max
	   if null, output results on cout
	*/
	public void CalcMinMaxAngle(double badellimit, Nullable<double>[] retvalues = null)
	{
	  int i;
	  int j;
	  int lpi1;
	  int lpi2;
	  int lpi3;
	  int lpi4;
	  double phimax = 0;
	  double phimin = 10;
	  double facephimax = 0;
	  double facephimin = 10;
	  int illegaltets = 0;
	  int negativetets = 0;
	  int badtets = 0;

	  for (i = 1; i <= GetNE(); i++)
	  {
		  int badel = 0;

		  Element el = VolumeElement(i);

		  if (el.GetType() != ELEMENT_TYPE.TET)
		  {
			  VolumeElement(i).flags.badel = 0;
			  continue;
		  }

		  if (el.Volume(Points()) < 0)
		  {
			  badel = 1;
			  negativetets++;
		  }


		  if (!LegalTet(el))
		  {
			  badel = 1;
			  illegaltets++;
			  (*testout) << "illegal tet: " << i << " ";
			  for (j = 1; j <= el.GetNP(); j++)
			  {
				(*testout) << el.PNum(j) << " ";
			  }
			  (*testout) << "\n";
		  }


		  // angles between faces
		  for (lpi1 = 1; lpi1 <= 3; lpi1++)
		  {
			for (lpi2 = lpi1 + 1; lpi2 <= 4; lpi2++)
			{
				lpi3 = 1;
				while (lpi3 == lpi1 || lpi3 == lpi2)
				{
				  lpi3++;
				}
				lpi4 = 10 - lpi1 - lpi2 - lpi3;

				Point3d p1 = new Point(el.PNum(lpi1));
				Point3d p2 = new Point(el.PNum(lpi2));
				Point3d p3 = new Point(el.PNum(lpi3));
				Point3d p4 = new Point(el.PNum(lpi4));

				Vec3d n = new Vec3d(p1, p2);
				n /= n.Length();
				Vec3d v1 = new Vec3d(p1, p3);
				Vec3d v2 = new Vec3d(p1, p4);

				v1 -= (n * v1) * n;
				v2 -= (n * v2) * n;

				double cosphi = (v1 * v2) / (v1.Length() * v2.Length());
				double phi = Math.Acos(cosphi);
				if (phi > phimax)
				{
					phimax = phi;
				}
				if (phi < phimin)
				{
					phimin = phi;
				}

				if ((180 / DefineConstants.M_PI) * phi > badellimit)
				{
				  badel = 1;
				}
			}
		  }


		  // angles in faces
		  for (j = 1; j <= 4; j++)
		  {
			  Element2d face = new Element2d(ELEMENT_TYPE.TRIG);
			  el.GetFace(j, face);
			  for (lpi1 = 1; lpi1 <= 3; lpi1++)
			  {
				  lpi2 = lpi1 % 3 + 1;
				  lpi3 = lpi2 % 3 + 1;

				  Point3d p1 = new Point(el.PNum(lpi1));
				  Point3d p2 = new Point(el.PNum(lpi2));
				  Point3d p3 = new Point(el.PNum(lpi3));

				  Vec3d v1 = new Vec3d(p1, p2);
				  Vec3d v2 = new Vec3d(p1, p3);
				  double cosphi = (v1 * v2) / (v1.Length() * v2.Length());
				  double phi = Math.Acos(cosphi);
				  if (phi > facephimax)
				  {
					  facephimax = phi;
				  }
				  if (phi < facephimin)
				  {
					  facephimin = phi;
				  }

				  if ((180 / DefineConstants.M_PI) * phi > badellimit)
				  {
					badel = 1;
				  }

			  }
		  }


		  VolumeElement(i).flags.badel = badel;
		  if (badel != 0)
		  {
			  badtets++;
		  }
	  }

	  if (GetNE() == null)
	  {
		  phimin = phimax = facephimin = facephimax = 0;
	  }

	  if (retvalues == null)
	  {
		  PrintMessage(1, "");
		  PrintMessage(1, "between planes:  phimin = ", (180 / DefineConstants.M_PI) * phimin, " phimax = ", (180 / DefineConstants.M_PI) * phimax);
		  PrintMessage(1, "inside planes:   phimin = ", (180 / DefineConstants.M_PI) * facephimin, " phimax = ", (180 / DefineConstants.M_PI) * facephimax);
		  PrintMessage(1, "");
	  }
	  else
	  {
		  retvalues[0] = (180 / DefineConstants.M_PI) * facephimin;
		  retvalues[1] = (180 / DefineConstants.M_PI) * facephimax;
		  retvalues[2] = (180 / DefineConstants.M_PI) * phimin;
		  retvalues[3] = (180 / DefineConstants.M_PI) * phimax;
	  }
	  PrintMessage(3, "negative tets: ", negativetets);
	  PrintMessage(3, "illegal tets:  ", illegaltets);
	  PrintMessage(3, "bad tets:      ", badtets);
	}

	/*
	  Marks elements which are dangerous to refine
	  return: number of illegal elements
	*/
	public int MarkIllegalElements()
	{
	  int cnt = 0;
	  foreach (var el in VolumeElements())
	  {
		if (!LegalTet(el))
		{
		  cnt++;
		}
	  }
	  return cnt;
	}

	/// orient surface mesh, for one sub-domain only
	public void SurfaceMeshOrientation()
	{
	  int i;
	  int j;
	  int nse = GetNSE();

	  BitArray used = new BitArray(nse);
	  used.Clear();
	  INDEX_2_HASHTABLE<int> edges = new INDEX_2_HASHTABLE<int>(nse+1);

	  bool haschanged = false;


	  Element2d tri = SurfaceElement(1);
	  for (j = 1; j <= 3; j++)
	  {
		  INDEX_2 i2 = new INDEX_2(tri.PNumMod(j), tri.PNumMod(j + 1));
		  edges.Set(i2, 1);
	  }
	  used.Set(1);

	  bool unused;
	  do
	  {
		  bool changed;
		  do
		  {
			  changed = false;
			  for (i = 1; i <= nse; i++)
			  {
				if (!used.Test(i))
				{
					Element2d el = surfelements.Elem(i);
					int found = 0;
					int foundrev = 0;
					for (j = 1; j <= 3; j++)
					{
						INDEX_2 i2 = new INDEX_2(el.PNumMod(j), el.PNumMod(j + 1));
						if (edges.Used(i2))
						{
						  foundrev = 1;
						}
						swap(i2.I1(), i2.I2());
						if (edges.Used(i2))
						{
						  found = 1;
						}
					}

					if (found != 0 || foundrev != 0)
					{
						if (foundrev != 0)
						{
						  swap(el.PNum(2), el.PNum(3));
						}

						changed = true;
						for (j = 1; j <= 3; j++)
						{
							INDEX_2 i2 = new INDEX_2(el.PNumMod(j), el.PNumMod(j + 1));
							edges.Set(i2, 1);
						}
						used.Set(i);
					}
				}
			  }
			  if (changed)
			  {
				haschanged = true;
			  }
		  } while (changed);


		  unused = false;
		  for (i = 1; i <= nse; i++)
		  {
			if (!used.Test(i))
			{
				unused = true;
				Element2d tri = SurfaceElement(i);
				for (j = 1; j <= 3; j++)
				{
					INDEX_2 i2 = new INDEX_2(tri.PNumMod(j), tri.PNumMod(j + 1));
					edges.Set(i2, 1);
				}
				used.Set(i);
				break;
			}
		  }
	  } while (unused);

	  if (haschanged)
	  {
		timestamp = netgen.GlobalMembers.NextTimeStamp();
	  }
	}

	/// convert mixed element mesh to tet-mesh
	public void Split2Tets()
	{
	  PrintMessage(1, "Split To Tets");
	  bool has_prisms = false;

	  int oldne = GetNE();
	  for (int i = 1; i <= oldne; i++)
	  {
		  Element el = VolumeElement(i);

		  if (el.GetType() == ELEMENT_TYPE.PRISM)
		  {
			  // prism, to 3 tets

			  // make minimal node to node 1
			  int minpi = 0;
			  PointIndex minpnum = new PointIndex();
			  minpnum = GetNP() + 1;

			  for (int j = 1; j <= 6; j++)
			  {
				  if (el.PNum(j) < minpnum)
				  {
					  minpnum.CopyFrom(el.PNum(j));
					  minpi = j;
				  }
			  }

			  if (minpi >= 4)
			  {
				  for (int j = 1; j <= 3; j++)
				  {
					swap(el.PNum(j), el.PNum(j + 3));
				  }
				  minpi -= 3;
			  }

			  while (minpi > 1)
			  {
				  int hi = 0;
				  for (int j = 0; j <= 3; j += 3)
				  {
					  hi = el.PNum(1 + j);
					  el.PNum(1 + j) = el.PNum(2 + j);
					  el.PNum(2 + j) = el.PNum(3 + j);
					  el.PNum(3 + j) = hi;
				  }
				  minpi--;
			  }

			  /*
			    version 1: edge from pi2 to pi6,
			    version 2: edge from pi3 to pi5,
			  */

			  int[][] ntets =
			  {
				  new int[] {1, 4, 5, 6, 1, 2, 3, 6, 1, 2, 5, 6},
				  new int[] {1, 4, 5, 6, 1, 2, 3, 5, 3, 1, 5, 6}
			  };

//C++ TO C# CONVERTER TODO TASK: C# does not have an equivalent to pointers to value types:
//ORIGINAL LINE: const int * min2pi;
			  int min2pi;

			  if (netgen.GlobalMembers.min2(el.PNum(2), el.PNum(6)) < netgen.GlobalMembers.min2(el.PNum(3), el.PNum(5)))
			  {
				  min2pi = ntets[0][0];
				  // (*testout) << "version 1 ";
			  }
			  else
			  {
				  min2pi = ntets[1][0];
				  // (*testout) << "version 2 ";
			  }


			  int firsttet = 1;
			  for (int j = 1; j <= 3; j++)
			  {
				  Element nel = new Element(ELEMENT_TYPE.TET);
				  for (int k = 1; k <= 4; k++)
				  {
					nel.PNum(k) = el.PNum(min2pi[4 * j + k - 5]);
				  }
				  nel.SetIndex(el.GetIndex());

				  int legal = 1;
				  for (int k = 1; k <= 3; k++)
				  {
					for (int l = k + 1; l <= 4; l++)
					{
					  if (nel.PNum(k) == nel.PNum(l))
					  {
						legal = 0;
					  }
					}
				  }

				  // (*testout) << nel << " ";
				  if (legal != 0)
				  {
					  if (firsttet != 0)
					  {
						  VolumeElement(i) = nel;
						  firsttet = 0;
					  }
					  else
					  {
						  AddVolumeElement(nel);
					  }
				  }
			  }
			  if (firsttet != 0)
			  {
				  Console.Write("no legal");
			  }
			  (*testout) << "\n";
		  }



		  else if (el.GetType() == ELEMENT_TYPE.HEX)
		  {
			  // hex to A) 2 prisms or B) to 5 tets

			  // make minimal node to node 1
			  int minpi = 0;
			  PointIndex minpnum = new PointIndex();
			  minpnum = GetNP() + 1;

			  for (int j = 1; j <= 8; j++)
			  {
				  if (el.PNum(j) < minpnum)
				  {
					  minpnum.CopyFrom(el.PNum(j));
					  minpi = j;
				  }
			  }

			  if (minpi >= 5)
			  {
				  for (int j = 1; j <= 4; j++)
				  {
					swap(el.PNum(j), el.PNum(j + 4));
				  }
				  minpi -= 4;
			  }

			  while (minpi > 1)
			  {
				  int hi = 0;
				  for (int j = 0; j <= 4; j += 4)
				  {
					  hi = el.PNum(1 + j);
					  el.PNum(1 + j) = el.PNum(2 + j);
					  el.PNum(2 + j) = el.PNum(3 + j);
					  el.PNum(3 + j) = el.PNum(4 + j);
					  el.PNum(4 + j) = hi;
				  }
				  minpi--;
			  }



			  int[][] to_prisms =
			  {
				  new int[] {0, 1, 2, 4, 5, 6, 0, 2, 3, 4, 6, 7},
				  new int[] {0, 1, 5, 3, 2, 6, 0, 5, 4, 3, 6, 7},
				  new int[] {0, 7, 4, 1, 6, 5, 0, 3, 7, 1, 2, 6}
			  };

			  int[] min2pi = 0;
			  if (netgen.GlobalMembers.min2(new netgen.Element(el[4]), new netgen.Element(el[6])) < netgen.GlobalMembers.min2(new netgen.Element(el[5]), new netgen.Element(el[7])))
			  {
				min2pi = to_prisms[0][0];
			  }
			  else if (netgen.GlobalMembers.min2(new netgen.Element(el[3]), new netgen.Element(el[6])) < netgen.GlobalMembers.min2(new netgen.Element(el[2]), new netgen.Element(el[7])))
			  {
				min2pi = to_prisms[1][0];
			  }
			  else if (netgen.GlobalMembers.min2(new netgen.Element(el[1]), new netgen.Element(el[6])) < netgen.GlobalMembers.min2(new netgen.Element(el[2]), new netgen.Element(el[5])))
			  {
				min2pi = to_prisms[2][0];
			  }

			  if (min2pi)
			  {
				  has_prisms = true;
				  for (int j = 0; j < 2; j++)
				  {
					  Element nel = new Element(ELEMENT_TYPE.PRISM);
					  for (int k = 0; k < 6; k++)
					  {
						nel[k] = el[min2pi[6 * j + k]];
					  }
					  nel.SetIndex(el.GetIndex());

					  if (j == 0)
					  {
						VolumeElement(i) = nel;
					  }
					  else
					  {
						AddVolumeElement(nel);
					  }
				  }
			  }
			  else
			  {
				  // split to 5 tets

				  int[] to_tets = {1, 2, 0, 5, 3, 0, 2, 7, 4, 5, 7, 0, 6, 7, 5, 2, 0, 2, 7, 5};

				  for (int j = 0; j < 5; j++)
				  {
					  Element nel = new Element(ELEMENT_TYPE.TET);
					  for (int k = 0; k < 4; k++)
					  {
						nel[k] = el[to_tets[4 * j + k]];
					  }
					  nel.SetIndex(el.GetIndex());

					  if (j == 0)
					  {
						VolumeElement(i) = nel;
					  }
					  else
					  {
						AddVolumeElement(nel);
					  }
				  }

			  }
		  }





		  else if (el.GetType() == ELEMENT_TYPE.PYRAMID)
		  {
			  // pyramid, to 2 tets

			  // cout << "pyramid: " << el << endl;

			  int[][] ntets =
			  {
				  new int[] {1, 2, 3, 5, 1, 3, 4, 5},
				  new int[] {1, 2, 4, 5, 4, 2, 3, 5}
			  };

			  int[] min2pi;

			  if (netgen.GlobalMembers.min2(new netgen.Element(el[0]), new netgen.Element(el[2])) < netgen.GlobalMembers.min2(new netgen.Element(el[1]), new netgen.Element(el[3])))
			  {
				min2pi = ntets[0][0];
			  }
			  else
			  {
				min2pi = ntets[1][0];
			  }

			  bool firsttet = true;
			  for (int j = 0; j < 2; j++)
			  {
				  Element nel = new Element(ELEMENT_TYPE.TET);
				  for (int k = 0; k < 4; k++)
				  {
					nel[k] = el[min2pi[4 * j + k] - 1];
				  }
				  nel.SetIndex(el.GetIndex());

				  // cout << "pyramid-tet: " << nel << endl;

				  bool legal = true;
				  for (int k = 0; k < 3; k++)
				  {
					for (int l = k + 1; l < 4; l++)
					{
					  if (nel[k] == nel[l])
					  {
						legal = false;
					  }
					}
				  }

				  if (legal)
				  {
					  (*testout) << nel << " ";
					  if (firsttet)
					  {
						VolumeElement(i) = nel;
					  }
					  else
					  {
						AddVolumeElement(nel);
					  }

					  firsttet = false;
				  }
			  }
			  if (firsttet)
			  {
				  Console.Write("no legal");
			  }
			  (*testout) << "\n";
		  }
	  }


	  int oldnse = GetNSE();
	  for (int i = 1; i <= oldnse; i++)
	  {
		  Element2d el = SurfaceElement(i);
		  if (el.GetNP() == 4)
		  {
			  (*testout) << "split el: " << el << " to ";

			  int[][] ntris =
			  {
				  new int[] {1, 2, 3, 1, 3, 4},
				  new int[] {1, 2, 4, 4, 2, 3}
			  };

			  int[] min2pi;

			  if (netgen.GlobalMembers.min2(el.PNum(1), el.PNum(3)) < netgen.GlobalMembers.min2(el.PNum(2), el.PNum(4)))
			  {
				min2pi = ntris[0][0];
			  }
			  else
			  {
				min2pi = ntris[1][0];
			  }

			  for (int j = 0; j < 6; j++)
			  {
				(*testout) << min2pi[j] << " ";
			  }


			  int firsttri = 1;
			  for (int j = 1; j <= 2; j++)
			  {
				  Element2d nel = new Element2d(3);
				  for (int k = 1; k <= 3; k++)
				  {
					nel.PNum(k) = el.PNum(min2pi[3 * j + k - 4]);
				  }
				  nel.SetIndex(el.GetIndex());

				  int legal = 1;
				  for (int k = 1; k <= 2; k++)
				  {
					for (int l = k + 1; l <= 3; l++)
					{
					  if (nel.PNum(k) == nel.PNum(l))
					  {
						legal = 0;
					  }
					}
				  }

				  if (legal != 0)
				  {
					  (*testout) << nel << " ";
					  if (firsttri != 0)
					  {
						  SurfaceElement(i) = nel;
						  firsttri = 0;
					  }
					  else
					  {
						  AddSurfaceElement(nel);
					  }
				  }
			  }
			  (*testout) << "\n";

		  }
	  }


	  if (has_prisms)
	  {

		Split2Tets();
	  }

	  else
	  {
		  for (int i = 1; i <= GetNE(); i++)
		  {
			  Element el = VolumeElement(i);
			  Point3d p1 = new Point(el.PNum(1));
			  Point3d p2 = new Point(el.PNum(2));
			  Point3d p3 = new Point(el.PNum(3));
			  Point3d p4 = new Point(el.PNum(4));

			  double vol = (new Vec3d(p1, p2) * netgen.GlobalMembers.Cross(new Vec3d(p1, p3), new Vec3d(p1, p4)));
			  if (vol > 0)
			  {
				swap(el.PNum(3), el.PNum(4));
			  }
		  }



		  UpdateTopology();
		  timestamp = netgen.GlobalMembers.NextTimeStamp();
	  }

	  RebuildSurfaceElementLists();
	}


	/// build box-search tree
	public void BuildElementSearchTree()
	{
	  if (elementsearchtreets == netgen.GlobalMembers.GetTimeStamp())
	  {
		  return;
	  }

	  lock (buildsearchtree_mutex)
	  {
		if (elementsearchtreets != netgen.GlobalMembers.GetTimeStamp())
		{
			NgLock @lock = new NgLock(object);
			@lock.Lock();

			PrintMessage(4, "Rebuild element searchtree");

			elementsearchtree = null;
			elementsearchtree = null;

			int ne = (dimension == 2) ? GetNSE() : GetNE();
			if (dimension == 3 && GetNE() == null && GetNSE() != null)
			{
			  ne = GetNSE();
			}

			if (ne != 0)
			{
				if (dimension == 2 || (dimension == 3 && GetNE() == null))
				{
					Box < 3> box(Box < 3>.EMPTY_BOX);
					for (SurfaceElementIndex sei = 0; sei < ne; sei++)
					{
					  box.Add(points[surfelements[sei].PNums()]);
					}

					box.Increase(1.01 * box.Diam());
					elementsearchtree = new BoxTree < 3> (box);

					for (SurfaceElementIndex sei = 0; sei < ne; sei++)
					{
						box.Set(points[surfelements[sei].PNums()]);
						elementsearchtree.Insert(box, sei + 1);
					}
				}
				else
				{
					Box < 3> box(Box < 3>.EMPTY_BOX);
					for (ElementIndex ei = 0; ei < ne; ei++)
					{
					  box.Add(points[volelements[ei].PNums()]);
					}

					box.Increase(1.01 * box.Diam());
					elementsearchtree = new BoxTree < 3> (box);

					for (ElementIndex ei = 0; ei < ne; ei++)
					{
						box.Set(points[volelements[ei].PNums()]);
						elementsearchtree.Insert(box, ei + 1);
					}
				}

				elementsearchtreets = netgen.GlobalMembers.GetTimeStamp();
			}
		}
	  }
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void SetPointSearchStartElement(const int el) const
	public void SetPointSearchStartElement(int el)
	{
		ps_startelement = el;
	}

	/// gives element of point, barycentric coordinates
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetElementOfPoint(const netgen::Point<3> & p, double lami[3] lami, bool build_searchtree = 0, const int index = -1, const bool allowindex = true) const
	public int GetElementOfPoint(netgen.Point < 3> p, double lami[3] lami, bool build_searchtree = false, int index = -1, bool allowindex = true)
	{
	  if (index != -1)
	  {
		  Array<int> dummy = new Array<int>(1);
		  dummy[0] = index;
		  return GetElementOfPoint(p.functorMethod, lami, dummy, build_searchtree, allowindex);
	  }
	  else
	  {
		return GetElementOfPoint(p.functorMethod,lami,null,build_searchtree,allowindex);
	  }
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetElementOfPoint(const netgen::Point<3> & p, double lami[3] lami, const Array<int> * const indices, bool build_searchtree = 0, const bool allowindex = true) const
	public int GetElementOfPoint(netgen.Point < 3> p, double lami[3] lami, Array<int> indices, bool build_searchtree = false, bool allowindex = true)
	{
	  // const double pointtol = 1e-12;
	  // netgen::Point<3> pmin = p - Vec<3> (pointtol, pointtol, pointtol);
	  // netgen::Point<3> pmax = p + Vec<3> (pointtol, pointtol, pointtol);

	  if (dimension == 2 || (dimension == 3 && GetNE() == null && GetNSE() != null))
	  {
		  int ne;
		  int ps_startelement = 0; // disable global buffering

		  if (ps_startelement != 0 && ps_startelement <= GetNSE() && PointContainedIn2DElement(p.functorMethod, lami, ps_startelement))
		  {
			return ps_startelement;
		  }

		  Array<int> locels = new Array<int>();
		  if (elementsearchtree || build_searchtree)
		  {
			  // update if necessary:
//C++ TO C# CONVERTER TODO TASK: There is no equivalent to 'const_cast' in C#:
			  const_cast<Mesh&>(this).BuildElementSearchTree();
			  // double tol = elementsearchtree->Tolerance();
			  // netgen::Point<3> pmin = p - Vec<3> (tol, tol, tol);
			  // netgen::Point<3> pmax = p + Vec<3> (tol, tol, tol);
			  elementsearchtree.GetIntersecting(p.functorMethod, p.functorMethod, locels);
			  ne = locels.Size();
		  }
		  else
		  {
			ne = GetNSE();
		  }

		  for (int i = 1; i <= ne; i++)
		  {
			  int ii;

			  if (elementsearchtree)
			  {
				ii = locels.Get(i);
			  }
			  else
			  {
				ii = i;
			  }

			  if (ii == ps_startelement)
			  {
				  continue;
			  }

			  if (indices != null && indices.Size() > 0)
			  {
				  bool contained = indices.Contains(SurfaceElement(ii).GetIndex());
				  if ((allowindex && !contained) || (!allowindex && contained))
				  {
					  continue;
				  }
			  }

			  if (PointContainedIn2DElement(p.functorMethod, lami, ii))
			  {
				  return ii;
			  }

		  }
		  return 0;
	  }
	  else

	  {
		  int ps_startelement = 0; // disable global buffering
		  // int i, j;
		  int ne;

		  if (ps_startelement != 0 && PointContainedIn3DElement(p.functorMethod, lami, ps_startelement))
		  {
			return ps_startelement;
		  }

		  Array<int> locels = new Array<int>();
		  if (elementsearchtree || build_searchtree)
		  {
			  // update if necessary:
//C++ TO C# CONVERTER TODO TASK: There is no equivalent to 'const_cast' in C#:
			  const_cast<Mesh&>(this).BuildElementSearchTree();
			  // double tol = elementsearchtree->Tolerance();
			  // netgen::Point<3> pmin = p - Vec<3> (tol, tol, tol);
			  // netgen::Point<3> pmax = p + Vec<3> (tol, tol, tol);
			  elementsearchtree.GetIntersecting(p.functorMethod, p.functorMethod, locels);
			  ne = locels.Size();
		  }
		  else
		  {
			ne = GetNE();
		  }

		  for (int i = 1; i <= ne; i++)
		  {
			  int ii;

			  if (elementsearchtree)
			  {
				ii = locels.Get(i);
			  }
			  else
			  {
				ii = i;
			  }
			  if (ii == ps_startelement)
			  {
				  continue;
			  }

			  if (indices != null && indices.Size() > 0)
			  {
				  bool contained = indices.Contains(VolumeElement(ii).GetIndex());
				  if ((allowindex && !contained) || (!allowindex && contained))
				  {
					  continue;
				  }
			  }

			  if (PointContainedIn3DElement(p.functorMethod, lami, ii))
			  {
				  ps_startelement = ii;
				  return ii;
			  }
		  }

		  // Not found, try uncurved variant:
		  for (int i = 1; i <= ne; i++)
		  {
			  int ii;

			  if (elementsearchtree)
			  {
				ii = locels.Get(i);
			  }
			  else
			  {
				ii = i;
			  }

			  if (indices != null && indices.Size() > 0)
			  {
				  bool contained = indices.Contains(VolumeElement(ii).GetIndex());
				  if ((allowindex && !contained) || (!allowindex && contained))
				  {
					  continue;
				  }
			  }


			  if (PointContainedIn3DElementOld(p.functorMethod, lami, ii))
			  {
				  ps_startelement = ii;
				  (*testout) << "WARNING: found element of point " << p.functorMethod << " only for uncurved mesh" << "\n";
				  return ii;
			  }
		  }


		  return 0;
	  }
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetSurfaceElementOfPoint(const netgen::Point<3> & p, double lami[3] lami, bool build_searchtree = 0, const int index = -1, const bool allowindex = true) const
	public int GetSurfaceElementOfPoint(netgen.Point < 3> p, double lami[3] lami, bool build_searchtree = false, int index = -1, bool allowindex = true)
	{
	  if (index != -1)
	  {
		  Array<int> dummy = new Array<int>(1);
		  dummy[0] = index;
		  return GetSurfaceElementOfPoint(p.functorMethod, lami, dummy, build_searchtree, allowindex);
	  }
	  else
	  {
		return GetSurfaceElementOfPoint(p.functorMethod,lami,null,build_searchtree,allowindex);
	  }
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetSurfaceElementOfPoint(const netgen::Point<3> & p, double lami[3] lami, const Array<int> * const indices, bool build_searchtree = 0, const bool allowindex = true) const
	public int GetSurfaceElementOfPoint(netgen.Point < 3> p, double lami[3] lami, Array<int> indices, bool build_searchtree = false, bool allowindex = true)
	{
	  if (dimension == 2)
	  {
		  throw new Exception("GetSurfaceElementOfPoint not yet implemented for 2D meshes");
	  }
	  else
	  {
		  double[] vlam = new double[3];
		  int velement = GetElementOfPoint(p.functorMethod,vlam,null,build_searchtree,allowindex);

		  //(*testout) << "p " << p << endl;
		  //(*testout) << "velement " << velement << endl;

		  if (GetNE() == null && GetNSE() != null)
		  {
			  lami[0] = vlam[0];
			  lami[1] = vlam[1];
			  lami[2] = vlam[2];
			  return velement;
		  }

		  Array<int> faces = new Array<int>();
		  topology.GetElementFaces(velement, faces);

		  //(*testout) << "faces " << faces << endl;

		  for (int i = 0; i < faces.Size(); i++)
		  {
			faces[i] = topology.GetFace2SurfaceElement(faces[i]);
		  }

		  //(*testout) << "surfel " << faces << endl;

		  for (int i = 0; i < faces.Size(); i++)
		  {
			  if (faces[i] == 0)
			  {
				continue;
			  }

			  if (indices != null && indices.Size() != 0)
			  {
				  if (indices.Contains(SurfaceElement(faces[i]).GetIndex()) && PointContainedIn2DElement(p.functorMethod, lami, faces[i], true))
				  {
					return faces[i];
				  }
			  }
			  else
			  {
				  if (PointContainedIn2DElement(p.functorMethod, lami, faces[i], true))
				  {
					  //(*testout) << "found point " << p << " in sel " << faces[i]
					  //	       << ", lam " << lami[0] << ", " << lami[1] << ", " << lami[2] << endl;
					  return faces[i];
				  }
			  }
		  }

		  Array<int> faces2 = new Array<int>();
		  topology.GetElementFaces(velement, faces2);
		  /*
		  cout << "no matching surf element" << endl
		       << "p = " << p << endl
		       << "faces-orig = " << faces2 << endl
		       << "faces = " << faces << endl
		       << ", vol el = " << velement
		       << ", vlam = " << vlam[0] << "," << vlam[1] << "," << vlam[2] << endl;
		  */
	  }

	  return 0;
	}

	/// give list of vol elements which are int the box(p1,p2)
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetIntersectingVolEls(const Point3d& p1, const Point3d& p2, Array<int> & locels) const
	public void GetIntersectingVolEls(Point3d p1, Point3d p2, Array<int> locels)
	{
	  elementsearchtree.GetIntersecting(p1, p2, locels);
	}

	///
	public int AddFaceDescriptor(FaceDescriptor fd)
	{
		facedecoding.Append(fd);
		return facedecoding.Size();
	}

	public int AddEdgeDescriptor(EdgeDescriptor fd)
	{
		edgedecoding.Append(fd);
		return edgedecoding.Size() - 1;
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: auto GetCommunicator() const
	public ngcore.NgMPI_Comm GetCommunicator()
	{
		return new ngcore.NgMPI_Comm(this.comm);
	}
	public void SetCommunicator(NgMPI_Comm acomm)
	{
	  this.comm = acomm;
	}

	///
	public void SetMaterial(int domnr, string mat)
	{
	  if (domnr > materials.Size())
	  {
		  int olds = materials.Size();
		  materials.SetSize(domnr);
		  for (int i = olds; i < domnr - 1; i++)
		  {
			materials[i] = "default";
		  }
	  }
	  /*
	  materials.Elem(domnr) = new char[strlen(mat)+1];
	  strcpy (materials.Elem(domnr), mat);
	  */
	  materials.Elem(domnr) = new string(mat);
	}

	///
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
//C++ TO C# CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
	string GetMaterial_emptystring("default");

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const string & GetMaterial(int domnr) const
	public string GetMaterial(int domnr)
	{
	  if (domnr <= materials.Size())
	  {
		return materials.Get(domnr);
	  }
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	  static string emptystring("default");
	  return GetMaterial_emptystring;
	}

	public static DLL_HEADER string defaultmat = new DLL_HEADER();
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const string * GetMaterialPtr(int domnr) const
	public string GetMaterialPtr(int domnr) // 1-based
	{
	  return domnr <= materials.Size() ? materials.Get(domnr) : GlobalMembers.defaultmat;
	}

	public void SetNBCNames(int nbcn)
	{
	  if (bcnames.Size())
	  {
		for (int i = 0; i < bcnames.Size(); i++)
		{
		  if (bcnames[i])
		  {
			  bcnames[i] = null;
		  }
		}
	  }
	  bcnames.SetSize(nbcn);
	  bcnames = 0;
	}

	public void SetBCName(int bcnr, string abcname)
	{
	  if (bcnr >= bcnames.Size())
	  {
		  int oldsize = bcnames.Size();
		  bcnames.SetSize(bcnr + 1); // keeps contents
		  for (int i = oldsize; i <= bcnr; i++)
		  {
			bcnames[i] = null;
		  }
	  }

	  if (bcnames[bcnr])
	  {
		  bcnames[bcnr] = null;
	  }
	  if (abcname != "default")
	  {
		bcnames[bcnr] = new string(abcname);
	  }
	  else
	  {
		bcnames[bcnr] = null;
	  }

	  foreach (var fd in facedecoding)
	  {
		if (fd.BCProperty() <= bcnames.Size())
		{
		  fd.SetBCName(bcnames[fd.BCProperty() - 1]);
		}
	  }
	}

//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	private string GetBCName_defaultstring = "default";

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const string & GetBCName(int bcnr) const
	public string GetBCName(int bcnr)
	{
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	  static string defaultstring = "default";

	  if (!bcnames.Size())
	  {
		return GetBCName_defaultstring;
	  }

	  if (bcnr < 0 || bcnr >= bcnames.Size())
	  {
		throw new Exception("illegal bc-number");
	  }

	  if (bcnames[bcnr])
	  {
		return bcnames[bcnr];
	  }
	  else
	  {
		return GetBCName_defaultstring;
	  }
	}

	public void SetNCD2Names(int ncd2n)
	{
	  if (cd2names.Size())
	  {
		for (int i = 0; i < cd2names.Size(); i++)
		{
	  if (cd2names[i])
	  {
		  cd2names[i] = null;
	  }
		}
	  }
	  cd2names.SetSize(ncd2n);
	  cd2names = 0;
	}

	public void SetCD2Name(int cd2nr, string abcname)
	{
	  cd2nr--;
	  (*testout) << "setCD2Name on edge " << cd2nr << " to " << abcname << "\n";
	  if (cd2nr >= cd2names.Size())
	  {
	  int oldsize = cd2names.Size();
	  cd2names.SetSize(cd2nr + 1);
	  for (int i = oldsize; i <= cd2nr; i++)
	  {
		cd2names[i] = null;
	  }
	  }
	  //if (cd2names[cd2nr]) delete cd2names[cd2nr];
	  if (abcname != "default")
	  {
		cd2names[cd2nr] = new string(abcname);
	  }
	  else
	  {
		cd2names[cd2nr] = null;
	  }
	}

//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	private string GetCD2Name_defaultstring = "default";

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const string & GetCD2Name(int cd2nr) const
	public string GetCD2Name(int cd2nr)
	{
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	  static string defaultstring = "default";
	  if (!cd2names.Size())
	  {
		return GetCD2Name_defaultstring;
	  }

	  if (cd2nr < 0 || cd2nr >= cd2names.Size())
	  {
		return GetCD2Name_defaultstring;
	  }

	  if (cd2names[cd2nr])
	  {
		return cd2names[cd2nr];
	  }
	  else
	  {
		return GetCD2Name_defaultstring;
	  }
	}

	public static DLL_HEADER string cd2_default_name = new DLL_HEADER();
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: string * GetCD2NamePtr(int cd2nr) const
	public string GetCD2NamePtr(int cd2nr)
	{
	  if (cd2nr < cd2names.Size() && cd2names[cd2nr])
	  {
		  return cd2names[cd2nr];
	  }
	  return GlobalMembers.cd2_default_name;
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: uint GetNCD2Names() const
	public uint GetNCD2Names()
	{
		return cd2names.Size();
	}

	public void SetNCD3Names(int ncd3n)
	{
	  if (cd3names.Size())
	  {
		for (int i = 0; i < cd3names.Size(); i++)
		{
	  if (cd3names[i])
	  {
		  cd3names[i] = null;
	  }
		}
	  }
	  cd3names.SetSize(ncd3n);
	  cd3names = 0;
	}

	public void SetCD3Name(int cd3nr, string abcname)
	{
	  cd3nr--;
	  (*testout) << "setCD3Name on vertex " << cd3nr << " to " << abcname << "\n";
	  if (cd3nr >= cd3names.Size())
	  {
	  int oldsize = cd3names.Size();
	  cd3names.SetSize(cd3nr + 1);
	  for (int i = oldsize; i <= cd3nr; i++)
	  {
		cd3names[i] = null;
	  }
	  }
	  if (abcname != "default")
	  {
		cd3names[cd3nr] = new string(abcname);
	  }
	  else
	  {
		cd3names[cd3nr] = null;
	  }
	}

//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	private string GetCD3Name_defaultstring = "default";

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const string & GetCD3Name(int cd3nr) const
	public string GetCD3Name(int cd3nr)
	{
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	  static string defaultstring = "default";
	  if (!cd3names.Size())
	  {
		return GetCD3Name_defaultstring;
	  }

	  if (cd3nr < 0 || cd3nr >= cd3names.Size())
	  {
		return GetCD3Name_defaultstring;
	  }

	  if (cd3names[cd3nr])
	  {
		return cd3names[cd3nr];
	  }
	  else
	  {
		return GetCD3Name_defaultstring;
	  }
	}

	public static DLL_HEADER string cd3_default_name = new DLL_HEADER();
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: string * GetCD3NamePtr(int cd3nr) const
	public string GetCD3NamePtr(int cd3nr)
	{
	  if (cd3nr < cd3names.Size() && cd3names[cd3nr])
	  {
		  return cd3names[cd3nr];
	  }
	  return GlobalMembers.cd3_default_name;
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: uint GetNCD3Names() const
	public uint GetNCD3Names()
	{
		return cd3names.Size();
	}

	public static DLL_HEADER string default_bc = new DLL_HEADER();
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: string * GetBCNamePtr(int bcnr) const
	public string GetBCNamePtr(int bcnr)
	{
		return (bcnr < bcnames.Size() && bcnames[bcnr]) ? bcnames[bcnr] : GlobalMembers.default_bc;
	}

	///
	public void ClearFaceDescriptors()
	{
		facedecoding.SetSize(0);
	}

	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetNFD() const
	public int GetNFD()
	{
		return facedecoding.Size();
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const FaceDescriptor & GetFaceDescriptor(int i) const
	public FaceDescriptor GetFaceDescriptor(int i)
	{
		return facedecoding.Get(i);
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const EdgeDescriptor & GetEdgeDescriptor(int i) const
	public EdgeDescriptor GetEdgeDescriptor(int i)
	{
		return edgedecoding[i];
	}


	///
	public FaceDescriptor GetFaceDescriptor(int i)
	{
		return facedecoding.Elem(i);
	}

	// #ifdef NONE
	//   /*
	//     Identify points pi1 and pi2, due to
	//     identification nr identnr
	//   */
	//   void AddIdentification (int pi1, int pi2, int identnr);

	//   int GetIdentification (int pi1, int pi2) const;
	//   int GetIdentificationSym (int pi1, int pi2) const;
	//   ///
	//   INDEX_2_HASHTABLE<int> & GetIdentifiedPoints () 
	//   { 
	//     return *identifiedpoints; 
	//   }

	//   ///
	//   void GetIdentificationMap (int identnr, Array<int> & identmap) const;
	//   ///
	//   void GetIdentificationPairs (int identnr, Array<INDEX_2> & identpairs) const;
	//   ///
	//   int GetMaxIdentificationNr () const
	//   { 
	//     return maxidentnr; 
	//   }
	// #endif

	/// return periodic, close surface etc. identifications
	public Identifications GetIdentifications()
	{
		return ident;
	}
	/// return periodic, close surface etc. identifications
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const Identifications & GetIdentifications() const
	public Identifications GetIdentifications()
	{
		return ident;
	}
	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool HasIdentifications() const
	public bool HasIdentifications()
	{
		return ident != null;
	}


	// #ifdef NONE
	//   void Mesh :: AddIdentification (int pi1, int pi2, int identnr)
	//   {
	//     INDEX_2 pair(pi1, pi2);
	//     //  pair.Sort();
	//     identifiedpoints->Set (pair, identnr);
	//     if (identnr > maxidentnr)
	//       maxidentnr = identnr;
	//     timestamp = NextTimeStamp();
	//   }

	//   int Mesh :: GetIdentification (int pi1, int pi2) const
	//   {
	//     INDEX_2 pair(pi1, pi2);
	//     if (identifiedpoints->Used (pair))
	//       return identifiedpoints->Get(pair);
	//     else
	//       return 0;
	//   }

	//   int Mesh :: GetIdentificationSym (int pi1, int pi2) const
	//   {
	//     INDEX_2 pair(pi1, pi2);
	//     if (identifiedpoints->Used (pair))
	//       return identifiedpoints->Get(pair);

	//     pair = INDEX_2 (pi2, pi1);
	//     if (identifiedpoints->Used (pair))
	//       return identifiedpoints->Get(pair);

	//     return 0;
	//   }


	//   void Mesh :: GetIdentificationMap (int identnr, Array<int> & identmap) const
	//   {
	//     int i, j;

	//     identmap.SetSize (GetNP());
	//     for (i = 1; i <= identmap.Size(); i++)
	//       identmap.Elem(i) = 0;

	//     for (i = 1; i <= identifiedpoints->GetNBags(); i++)
	//       for (j = 1; j <= identifiedpoints->GetBagSize(i); j++)
	// 	{
	// 	  INDEX_2 i2;
	// 	  int nr;
	// 	  identifiedpoints->GetData (i, j, i2, nr);

	// 	  if (nr == identnr)
	// 	    {
	// 	      identmap.Elem(i2.I1()) = i2.I2();
	// 	    }
	// 	}
	//   }


	//   void Mesh :: GetIdentificationPairs (int identnr, Array<INDEX_2> & identpairs) const
	//   {
	//     int i, j;

	//     identpairs.SetSize(0);

	//     for (i = 1; i <= identifiedpoints->GetNBags(); i++)
	//       for (j = 1; j <= identifiedpoints->GetBagSize(i); j++)
	// 	{
	// 	  INDEX_2 i2;
	// 	  int nr;
	// 	  identifiedpoints->GetData (i, j, i2, nr);

	// 	  if (identnr == 0 || nr == identnr)
	// 	    identpairs.Append (i2);
	// 	}
	//   }
	// #endif



//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void InitPointCurve(double red = 1, double green = 0, double blue = 0) const
	public void InitPointCurve(double red = 1, double green = 0, double blue = 0)
	{
	  pointcurves_startpoint.Append(pointcurves.Size());
	  pointcurves_red.Append(red);
	  pointcurves_green.Append(green);
	  pointcurves_blue.Append(blue);
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void AddPointCurvePoint(const Point3d & pt) const
	public void AddPointCurvePoint(Point3d pt)
	{
	  pointcurves.Append(pt);
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetNumPointCurves() const
	public int GetNumPointCurves()
	{
	  return pointcurves_startpoint.Size();
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetNumPointsOfPointCurve(int curve) const
	public int GetNumPointsOfPointCurve(int curve)
	{
	  if (curve == pointcurves_startpoint.Size() - 1)
	  {
		return (pointcurves.Size() - pointcurves_startpoint.Last());
	  }
	  else
	  {
		return (pointcurves_startpoint[curve+1] - pointcurves_startpoint[curve]);
	  }
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: Point3d & GetPointCurvePoint(int curve, int n) const
	public Point3d GetPointCurvePoint(int curve, int n)
	{
	  return pointcurves[pointcurves_startpoint[curve] + n];
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetPointCurveColor(int curve, double & red, double & green, double & blue) const
	public void GetPointCurveColor(int curve, ref double red, ref double green, ref double blue)
	{
	  red = pointcurves_red[curve];
	  green = pointcurves_green[curve];
	  blue = pointcurves_blue[curve];
	}




	/// find number of vertices
	public void ComputeNVertices()
	{
	  numvertices = 0;

	  foreach (Element el in VolumeElements())
	  {
		foreach (PointIndex v in el.Vertices())
		{
		  if (v > numvertices)
		  {
			  numvertices = v;
		  }
		}
	  }

	  foreach (Element2d el in SurfaceElements())
	  {
		foreach (PointIndex v in el.Vertices())
		{
		  if (v > numvertices)
		  {
			  numvertices = v;
		  }
		}
	  }

	  numvertices += 1 - PointIndex.BASE;
	}

	/// number of vertices (no edge-midpoints)
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetNV() const
	public int GetNV()
	{
	  if (numvertices < 0)
	  {
		return new auto(GetNP());
	  }
	  else
	  {
		return numvertices;
	  }
	}

	/// remove edge points
	public void SetNP(int np)
	{
	  points.SetSize(np);
	  //  ptyps.SetSize(np);

	  int mlold = mlbetweennodes.Size();
	  mlbetweennodes.SetSize(np);
	  if (np > mlold)
	  {
		for (int i = mlold + PointIndex.BASE; i < np + PointIndex.BASE; i++)
		{
			mlbetweennodes[i].I1() = PointIndex.BASE-1;
			mlbetweennodes[i].I2() = PointIndex.BASE-1;
		}
	  }

	  GetIdentifications().SetMaxPointNr(np + PointIndex.BASE-1);
	}





	/*
	  void Mesh :: BuildConnectedNodes ()
	  {
	  if (PureTetMesh())
	  {
	  connectedtonode.SetSize(0);
	  return;
	  }
  
  
	  int i, j, k;
	  int np = GetNP();
	  int ne = GetNE();
	  TABLE<int> conto(np);
	  for (i = 1; i <= ne; i++)
	  {
	  const Element & el = VolumeElement(i);
  
	  if (el.GetType() == PRISM)
	  {
	  for (j = 1; j <= 6; j++)
	  {
	  int n1 = el.PNum (j);
	  int n2 = el.PNum ((j+2)%6+1);
	  //	    if (n1 != n2)
	  {
	  int found = 0;
	  for (k = 1; k <= conto.EntrySize(n1); k++)
	  if (conto.Get(n1, k) == n2)
	  {
	  found = 1;
	  break;
	  }
	  if (!found)
	  conto.Add (n1, n2);
	  }
	  }
	  }
	  else if (el.GetType() == PYRAMID)
	  {
	  for (j = 1; j <= 4; j++)
	  {
	  int n1, n2;
	  switch (j)
	  {
	  case 1: n1 = 1; n2 = 4; break;
	  case 2: n1 = 4; n2 = 1; break;
	  case 3: n1 = 2; n2 = 3; break;
	  case 4: n1 = 3; n2 = 2; break;
	  }
  
	  int found = 0;
	  for (k = 1; k <= conto.EntrySize(n1); k++)
	  if (conto.Get(n1, k) == n2)
	  {
	  found = 1;
	  break;
	  }
	  if (!found)
	  conto.Add (n1, n2);
	  }
	  }
	  }
  
	  connectedtonode.SetSize(np);
	  for (i = 1; i <= np; i++)
	  connectedtonode.Elem(i) = 0;
  
	  for (i = 1; i <= np; i++)
	  if (connectedtonode.Elem(i) == 0)
	  {
	  connectedtonode.Elem(i) = i;
	  ConnectToNodeRec (i, i, conto);
	  }
  
  
  
	  }
  
	  void Mesh :: ConnectToNodeRec (int node, int tonode, 
	  const TABLE<int> & conto)
	  {
	  int i, n2;
	  //  (*testout) << "connect " << node << " to " << tonode << endl;
	  for (i = 1; i <= conto.EntrySize(node); i++)
	  {
	  n2 = conto.Get(node, i);
	  if (!connectedtonode.Get(n2))
	  {
	  connectedtonode.Elem(n2) = tonode;
	  ConnectToNodeRec (n2, tonode, conto);
	  }
	  }
	  }
	*/


//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool PureTrigMesh(int faceindex = 0) const
	public bool PureTrigMesh(int faceindex = 0)
	{
	  // if (!faceindex) return !mparam.quad;

	  if (faceindex == 0)
	  {
	  for (int i = 1; i <= GetNSE(); i++)
	  {
		if (SurfaceElement(i).GetNP() != 3)
		{
		  return false;
		}
	  }
	  return true;
	  }

	  for (int i = 1; i <= GetNSE(); i++)
	  {
		if (SurfaceElement(i).GetIndex() == faceindex && SurfaceElement(i).GetNP() != 3)
		{
		  return false;
		}
	  }
	  return true;
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool PureTetMesh() const
	public bool PureTetMesh()
	{
	  for (ElementIndex ei = 0; ei < GetNE(); ei++)
	  {
//C++ TO C# CONVERTER TODO TASK: The following line was determined to contain a copy constructor call - this should be verified and a copy constructor should be created:
//ORIGINAL LINE: if (VolumeElement(ei).GetNP() != 4)
		if (VolumeElement(new netgen.ElementIndex(ei)).GetNP() != 4)
		{
		  return false;
		}
	  }
	  return true;
	}


//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const MeshTopology & GetTopology() const
	public MeshTopology GetTopology()
	{
		return new netgen.MeshTopology(topology);
	}

	public void UpdateTopology(TaskManager tm = DummyTaskManager, Tracer tracer = DummyTracer)
	{
	  topology.Update(tm, tracer);
	  tracer("call update clusters", false);
	  clusters.Update(tm, tracer);
	  tracer("call update clusters", true);
#if PARALLEL
	  if (paralleltop != null)
	  {
		  paralleltop.Reset();
		  paralleltop.UpdateCoarseGrid();
	  }
#endif
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: class CurvedElements & GetCurvedElements() const
	public CurvedElements GetCurvedElements()
	{
		return curvedelems;
	}

//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//	DLL_HEADER void BuildCurvedElements(Refinement @ref, int aorder, bool arational = false);
	public void BuildCurvedElements(int aorder)
	{
	  if (GetGeometry() == null)
	  {
		throw new Exception("don't have a geometry for mesh curving");
	  }

	  GetCurvedElements().BuildCurvedElements(GetGeometry().GetRefinement(), aorder, false);

	  for (SegmentIndex seg = 0; seg < GetNSeg(); seg++)
	  {
//C++ TO C# CONVERTER TODO TASK: The following line was determined to contain a copy constructor call - this should be verified and a copy constructor should be created:
//ORIGINAL LINE: (*this)[seg].SetCurved(GetCurvedElements().IsSegmentCurved(seg));
		this[seg].SetCurved(GetCurvedElements().IsSegmentCurved(new netgen.SegmentIndex(seg)));
	  }
	  for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
	  {
//C++ TO C# CONVERTER TODO TASK: The following line was determined to contain a copy constructor call - this should be verified and a copy constructor should be created:
//ORIGINAL LINE: (*this)[sei].SetCurved(GetCurvedElements().IsSurfaceElementCurved(sei));
		this[sei].SetCurved(GetCurvedElements().IsSurfaceElementCurved(new netgen.SurfaceElementIndex(sei)));
	  }
	  for (ElementIndex ei = 0; ei < GetNE(); ei++)
	  {
//C++ TO C# CONVERTER TODO TASK: The following line was determined to contain a copy constructor call - this should be verified and a copy constructor should be created:
//ORIGINAL LINE: (*this)[ei].SetCurved(GetCurvedElements().IsElementCurved(ei));
		this[ei].SetCurved(GetCurvedElements().IsElementCurved(new netgen.ElementIndex(ei)));
	  }

	  SetNextMajorTimeStamp();
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const class AnisotropicClusters & GetClusters() const
	public AnisotropicClusters GetClusters()
	{
		return clusters;
	}


	public class CSurfaceArea
	{
	  private readonly Mesh mesh;
	  private bool valid;
	  private double area;
	  public CSurfaceArea(Mesh amesh)
	  {
		  this.mesh = new netgen.Mesh(amesh);
		  this.valid = false;
		  ;
	  }

	  public void Add(Element2d sel)
	  {
	if (sel.GetNP() == 3)
	{
	  area += netgen.GlobalMembers.Cross(mesh[sel[1]] - mesh[sel[0]], mesh[sel[2]] - mesh[sel[0]]).Length() / 2;
	}
	else
	{
	  area += netgen.GlobalMembers.Cross(new Vec3d(mesh[sel.PNum(1)], mesh[sel.PNum(3)]), new Vec3d(mesh[sel.PNum(1)], mesh[sel.PNum(4)])).Length() / 2;
	}
	  }
	  public void ReCalc()
	  {
	area = 0;
	for (SurfaceElementIndex sei = 0; sei < mesh.GetNSE(); sei++)
	{
	  Add(mesh[sei]);
	}
	valid = true;
	  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: operator double() const
	  public static implicit operator double(CSurfaceArea ImpliedObject)
	  {
		  return ImpliedObject.area;
	  }
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool Valid() const
	  public bool Valid()
	  {
		  return valid;
	  }
	}

	public CSurfaceArea surfarea = new CSurfaceArea();
	public CSurfaceArea SurfaceArea()
	{
//C++ TO C# CONVERTER TODO TASK: The following line was determined to contain a copy constructor call - this should be verified and a copy constructor should be created:
//ORIGINAL LINE: return surfarea;
		return new netgen.Mesh.CSurfaceArea(surfarea);
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const CSurfaceArea & SurfaceArea() const
	public CSurfaceArea SurfaceArea()
	{
//C++ TO C# CONVERTER TODO TASK: The following line was determined to contain a copy constructor call - this should be verified and a copy constructor should be created:
//ORIGINAL LINE: return surfarea;
		return new netgen.Mesh.CSurfaceArea(surfarea);
	}



//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetTimeStamp() const
	public int GetTimeStamp()
	{
		return timestamp;
	}
	public void SetNextTimeStamp()
	{
		timestamp = netgen.GlobalMembers.NextTimeStamp();
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetMajorTimeStamp() const
	public int GetMajorTimeStamp()
	{
		return majortimestamp;
	}
	public void SetNextMajorTimeStamp()
	{
		majortimestamp = timestamp = netgen.GlobalMembers.NextTimeStamp();
	}


	/// return mutex
	public NgMutex Mutex()
	{
		return object;
	}
	public NgMutex MajorMutex()
	{
//C++ TO C# CONVERTER TODO TASK: The following line was determined to contain a copy constructor call - this should be verified and a copy constructor should be created:
//ORIGINAL LINE: return majormutex;
		return new netgen.NgMutex(majormutex);
	}


//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: NetgenGeometry GetGeometry() const
	public NetgenGeometry GetGeometry()
	{
	  return geometry;
	}
	public void SetGeometry(NetgenGeometry geom)
	{
	  geometry = geom;
	}

	///
	public void SetUserData(string id, Array<int> data)
	{
	  if (userdata_int.Used(id))
	  {
		if (userdata_int[id] != null)
		{
			userdata_int[id].Dispose();
		}
	  }

	  Array<int> newdata = new Array<int>(data);

	  userdata_int.Set(id, newdata);
	}

	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool GetUserData(const char * id, Array<int> & data, int shift = 0) const
	public bool GetUserData(string id, Array<int> data, int shift = 0)
	{
	  if (userdata_int.Used(id))
	  {
		  if (data.Size() < (*userdata_int[id]).Size() + shift)
		  {
			data.SetSize((*userdata_int[id]).Size() + shift);
		  }
		  for (int i = 0; i < (*userdata_int[id]).Size(); i++)
		  {
			data[i + shift] = (*userdata_int[id])[i];
		  }
		  return true;
	  }
	  else
	  {
		  data.SetSize(0);
		  return false;
	  }
	}

	///
	public void SetUserData(string id, Array<double> data)
	{
	  if (userdata_double.Used(id))
	  {
		if (userdata_double[id] != null)
		{
			userdata_double[id].Dispose();
		}
	  }

	  Array<double> newdata = new Array<double>(data);

	  userdata_double.Set(id, newdata);
	}

	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool GetUserData(const char * id, Array<double> & data, int shift = 0) const
	public bool GetUserData(string id, Array<double> data, int shift = 0)
	{
	  if (userdata_double.Used(id))
	  {
		  if (data.Size() < (*userdata_double[id]).Size() + shift)
		  {
			data.SetSize((*userdata_double[id]).Size() + shift);
		  }
		  for (int i = 0; i < (*userdata_double[id]).Size(); i++)
		  {
			data[i + shift] = (*userdata_double[id])[i];
		  }
		  return true;
	  }
	  else
	  {
		  data.SetSize(0);
		  return false;
	  }
	}

	///
//C++ TO C# CONVERTER TODO TASK: C# has no concept of a 'friend' function:
//ORIGINAL LINE: friend void OptimizeRestart(Mesh & mesh3d);
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//	void OptimizeRestart(Mesh mesh3d);
	///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void PrintMemInfo(ostream & ost) const
	public void PrintMemInfo(ostream ost)
	{
	  ost << "Mesh Mem:" << "\n";

	  ost << GetNP() << " Points, of size " << sizeof(Point3d) << " + " << sizeof(POINTTYPE) << " = " << GetNP() * (sizeof(Point3d) + sizeof(POINTTYPE)) << "\n";

	  ost << GetNSE() << " Surface elements, of size " << sizeof(Element2d) << " = " << GetNSE() * sizeof(Element2d) << "\n";

	  ost << GetNE() << " Volume elements, of size " << sizeof(Element) << " = " << GetNE() * sizeof(Element) << "\n";

	  // ost << "surfs on node:";
	  // surfacesonnode.PrintMemInfo (cout);

	  ost << "boundaryedges: ";
	  if (boundaryedges != null)
	  {
		boundaryedges.PrintMemInfo(cout);
	  }

	  ost << "surfelementht: ";
	  if (surfelementht != null)
	  {
		surfelementht.PrintMemInfo(cout);
	  }
	}

	/// 
//C++ TO C# CONVERTER TODO TASK: C# has no concept of a 'friend' class:
//	friend class Meshing3;


	public enum GEOM_TYPE
	{
		NO_GEOM = 0,
		GEOM_2D = 1,
		GEOM_CSG = 10,
		GEOM_STL = 11,
		GEOM_OCC = 12,
		GEOM_ACIS = 13
	}
	public GEOM_TYPE geomtype;



#if PARALLEL
	/// returns parallel topology
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: class ParallelMeshTopology & GetParallelTopology() const
	public ParallelMeshTopology GetParallelTopology()
	{
		return paralleltop;
	}


	/// distributes the master-mesh to local meshes
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//	void Distribute();
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//	void Distribute(Array<int> volume_weights, Array<int> surface_weights, Array<int> segment_weights);


	/// find connection to parallel meshes
	//   void FindExchangePoints () ;

	//   void FindExchangeEdges ();
	//   void FindExchangeFaces ();

	/// use metis to decompose master mesh 
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//	void ParallelMetis(); //  Array<int> & neloc );
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//	void ParallelMetis(Array<int> volume_weights, Array<int> surface_weights, Array<int> segment_weights);

//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//	void PartHybridMesh(); //  Array<int> & neloc );
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//	void PartDualHybridMesh(); //  Array<int> & neloc );
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//	void PartDualHybridMesh2D(); // ( Array<int> & neloc );


	/// send mesh from master to local procs
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//	void SendRecvMesh();

	/// send mesh to parallel machine, keep global mesh at master 
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void SendMesh() const;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//	void SendMesh(); // Mesh * mastermesh, Array<int> & neloc) const;
	/// loads a mesh sent from master processor
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//	void ReceiveParallelMesh();

	public Array<int> vol_partition = new Array<int>();
	public Array<int> surf_partition = new Array<int>();
	public Array<int> seg_partition = new Array<int>();

#else
	public void Distribute()
	{
	}
	public void SendRecvMesh()
	{
	}
	public void Distribute(Array<int> volume_weights, Array<int> surface_weights, Array<int> segment_weights)
	{
	}
#endif


  }

}





