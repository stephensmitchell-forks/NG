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
/* File:   curvedelems.hpp                                                */
/* Author: Robert Gaisbauer (first version)                               */
/*         redesign by Joachim Schoeberl                                  */
/* Date:   27. Sep. 02, Feb 2006                                          */
/**************************************************************************/




//C++ TO C# CONVERTER NOTE: C# has no need of forward class declarations:
//class Refinement;


public class CurvedElements
{
  private readonly Mesh mesh;

  private Array<int> edgeorder = new Array<int>();
  private Array<int> faceorder = new Array<int>();

  private Array<int> edgecoeffsindex = new Array<int>();
  private Array<int> facecoeffsindex = new Array<int>();

  private Array< Vec < 3>> edgecoeffs = new Array< Vec < 3>>();
  private Array< Vec < 3>> facecoeffs = new Array< Vec < 3>>();

  private Array< double > edgeweight = new Array< double >(); // for rational 2nd order splines

  private int order;
  private bool rational;

  private bool ishighorder;

  public CurvedElements(Mesh amesh)
  {
	  this.mesh = amesh;
	order = 1;
	rational = 0;
	ishighorder = 0;
  }

  public void Dispose()
  {
	GlobalMembers.jacpols2.SetSize(0);
  }

  // bool IsHighOrder() const { return order > 1; }
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool IsHighOrder() const
  public bool IsHighOrder()
  {
	  return ishighorder;
  }

  // void SetHighOrder (int aorder) { order=aorder; }
  public void SetIsHighOrder(bool ho)
  {
	  ishighorder = ho;
  }

//C++ TO C# CONVERTER NOTE: Enums must be named in C#, so the following enum has been named by the converter:
  public enum AnonymousEnum
  {
	  MPI_TAG_CURVE = MPI_TAG_MESH + 20
  }
  public void BuildCurvedElements(Refinement @ref, int aorder, bool arational = false)
  {

	ishighorder = 0;
	order = 1;

	// MPI_Comm curve_comm;
	var curve_comm = mesh.GetCommunicator();
#if PARALLEL

	ParallelMeshTopology partop = mesh.GetParallelTopology();
	// MPI_Comm_dup (mesh.GetCommunicator(), &curve_comm);
	Array<int> procs = new Array<int>();
#else
	// curve_comm = mesh.GetCommunicator();
#endif
	int id = curve_comm.Rank();
	int ntasks = curve_comm.Size();

	bool working = (ntasks == 1) || (id > 0);

	if (working)
	{
	  order = aorder;
	}

	if (mesh.coarsemesh)
	{
	mesh.coarsemesh.GetCurvedElements().BuildCurvedElements(@ref, aorder, arational);
		order = aorder;
		rational = arational;
		ishighorder = (order > 1);
	return;
	}


	PrintMessage(1, "Curve elements, order = ", aorder);
	if (rational)
	{
		PrintMessage(1, "curved elements with rational splines");
	}

	// if (working)
//C++ TO C# CONVERTER TODO TASK: There is no equivalent to 'const_cast' in C#:
	const_cast<Mesh&> (mesh).UpdateTopology();
	MeshTopology top = mesh.GetTopology();

	rational = arational;

	Array<int> edgenrs = new Array<int>();
	int nedges = top.GetNEdges();
	int nfaces = top.GetNFaces();

	edgeorder.SetSize(nedges);
	faceorder.SetSize(nfaces);

	edgeorder = 1;
	faceorder = 1;

	if (rational)
	{
		edgeweight.SetSize(nedges);
		edgeweight = 1.0;
	}


	if (aorder <= 1)
	{
	for (ElementIndex ei = 0; ei < mesh.GetNE(); ei++)
	{
	  if (mesh[ei].GetType() == TET10)
	  {
		ishighorder = 1;
	  }
	}
	return;
	}


	if (rational)
	{
		aorder = 2;
	}

	if (working)
	{
	if (mesh.GetDimension() == 3)
	{
	  for (SurfaceElementIndex i = 0; i < mesh.GetNSE(); i++)
	  {
		  top.GetEdges(i, edgenrs);
		  for (int j = 0; j < edgenrs.Size(); j++)
		  {
		edgeorder[edgenrs[j]] = aorder;
		  }
		  faceorder[top.GetFace(i)] = aorder;
	  }
	}
	for (SegmentIndex i = 0; i < mesh.GetNSeg(); i++)
	{
	  edgeorder[top.GetEdge(i)] = aorder;
	}
	}

	if (rational)
	{
		edgeorder = 2;
		faceorder = 1;
	}


#if PARALLEL
	TABLE<int> send_orders = new TABLE<int>(ntasks);
	TABLE<int> recv_orders = new TABLE<int>(ntasks);

	if (ntasks > 1 && working)
	{
	for (int e = 0; e < edgeorder.Size(); e++)
	{
		partop.GetDistantEdgeNums(e+1, procs);
		for (int j = 0; j < procs.Size(); j++)
		{
		  send_orders.Add(procs[j], edgeorder[e]);
		}
	}
	for (int f = 0; f < faceorder.Size(); f++)
	{
		partop.GetDistantFaceNums(f + 1, procs);
		for (int j = 0; j < procs.Size(); j++)
		{
		  send_orders.Add(procs[j], faceorder[f]);
		}
	}
	}

	if (ntasks > 1)
	{
	  netgen.GlobalMembers.MyMPI_ExchangeTable(send_orders, recv_orders, MPI_TAG_CURVE, curve_comm);
	}

	if (ntasks > 1 && working)
	{
	Array<int> cnt = new Array<int>(ntasks);
	cnt = 0;
	for (int e = 0; e < edgeorder.Size(); e++)
	{
		partop.GetDistantEdgeNums(e+1, procs);
		for (int j = 0; j < procs.Size(); j++)
		{
		  edgeorder[e] = Math.Max(edgeorder[e], recv_orders[procs[j]][cnt[procs[j]]++]);
		}
	}
	for (int f = 0; f < faceorder.Size(); f++)
	{
		partop.GetDistantFaceNums(f + 1, procs);
		for (int j = 0; j < procs.Size(); j++)
		{
		  faceorder[f] = Math.Max(faceorder[f], recv_orders[procs[j]][cnt[procs[j]]++]);
		}
	}
	}
#endif


	edgecoeffsindex.SetSize(nedges + 1);
	int nd = 0;
	for (int i = 0; i < nedges; i++)
	{
	edgecoeffsindex[i] = nd;
	nd += Math.Max(0, edgeorder[i] - 1);
	}
	edgecoeffsindex[nedges] = nd;

	edgecoeffs.SetSize(nd);
	edgecoeffs = Vec < 3> (0,0,0);


	facecoeffsindex.SetSize(nfaces + 1);
	nd = 0;
	for (int i = 0; i < nfaces; i++)
	{
	facecoeffsindex[i] = nd;
	if (top.GetFaceType(i + 1) == TRIG)
	{
	  nd += netgen.GlobalMembers.max2(0, (faceorder[i] - 1) * (faceorder[i] - 2) / 2);
	}
	else
	{
	  nd += netgen.GlobalMembers.max2(0, netgen.GlobalMembers.sqr(faceorder[i] - 1));
	}
	}
	facecoeffsindex[nfaces] = nd;

	facecoeffs.SetSize(nd);
	facecoeffs = Vec < 3> (0,0,0);


	if (@ref == null || aorder <= 1)
	{
		order = aorder;
	return;
	}

	Array<double> xi = new Array<double>();
	Array<double> weight = new Array<double>();

	netgen.GlobalMembers.ComputeGaussRule(aorder + 4, xi, weight); // on (0,1)

	if (!GlobalMembers.jacpols2.Size())
	{
	GlobalMembers.jacpols2.SetSize(100);
	for (int i = 0; i < 100; i++)
	{
	  GlobalMembers.jacpols2[i] = new JacobiRecPol (100, i, 2);
	}
	}

	PrintMessage(3, "Curving edges");

	if (mesh.GetDimension() == 3 || rational)
	{
	Array<int> surfnr = new Array<int>(nedges);
	Array<PointGeomInfo> gi0 = new Array<PointGeomInfo>(nedges);
	Array<PointGeomInfo> gi1 = new Array<PointGeomInfo>(nedges);
	surfnr = -1;

	if (working)
	{
	  for (SurfaceElementIndex i = 0; i < mesh.GetNSE(); i++)
	  {
		  top.GetEdges(i, edgenrs);
		  Element2d el = mesh[i];
		  ELEMENT_EDGE edges = MeshTopology.GetEdges0(el.GetType());

		  for (int i2 = 0; i2 < edgenrs.Size(); i2++)
		  {
		  // PointIndex pi1 = el[edges[i2][0]];
		  // PointIndex pi2 = el[edges[i2][1]];

		  // bool swap = pi1 > pi2;

		  // Point<3> p1 = mesh[pi1];
		  // Point<3> p2 = mesh[pi2];

		  // int order1 = edgeorder[edgenrs[i2]];
		  // int ndof = max (0, order1-1);

		  surfnr[edgenrs[i2]] = mesh.GetFaceDescriptor(el.GetIndex()).SurfNr();
		  gi0[edgenrs[i2]] = el.GeomInfoPi(edges[i2][0] + 1);
		  gi1[edgenrs[i2]] = el.GeomInfoPi(edges[i2][1] + 1);
		  }
	  }
	}


#if PARALLEL
	if (ntasks > 1)
	{
		// distribute it ...
		TABLE<double> senddata = new TABLE<double>(ntasks);
		TABLE<double> recvdata = new TABLE<double>(ntasks);
		if (working)
		{
		  for (int e = 0; e < nedges; e++)
		  {
		  partop.GetDistantEdgeNums(e+1, procs);
		  for (int j = 0; j < procs.Size(); j++)
		  {
			  senddata.Add(procs[j], surfnr[e]);
			  if (surfnr[e] != -1)
			  {
			  senddata.Add(procs[j], gi0[e].trignum);
			  senddata.Add(procs[j], gi0[e].u);
			  senddata.Add(procs[j], gi0[e].v);
			  senddata.Add(procs[j], gi1[e].trignum);
			  senddata.Add(procs[j], gi1[e].u);
			  senddata.Add(procs[j], gi1[e].v);
			  }
		  }
		  }
		}

		netgen.GlobalMembers.MyMPI_ExchangeTable(senddata, recvdata, MPI_TAG_CURVE, curve_comm);


		Array<int> cnt = new Array<int>(ntasks);
		cnt = 0;
		if (working)
		{
		  for (int e = 0; e < nedges; e++)
		  {
		  partop.GetDistantEdgeNums(e+1, procs);
		  for (int j = 0; j < procs.Size(); j++)
		  {
			  int surfnr1 = recvdata[procs[j]][cnt[procs[j]]++];
			  if (surfnr1 != -1)
			  {
			  surfnr[e] = surfnr1;
			  gi0[e].trignum = (int)recvdata[procs[j]][cnt[procs[j]]++];
			  gi0[e].u = recvdata[procs[j]][cnt[procs[j]]++];
			  gi0[e].v = recvdata[procs[j]][cnt[procs[j]]++];
			  gi1[e].trignum = (int)recvdata[procs[j]][cnt[procs[j]]++];
			  gi1[e].u = recvdata[procs[j]][cnt[procs[j]]++];
			  gi1[e].v = recvdata[procs[j]][cnt[procs[j]]++];
			  }
		  }
		  }
		}

	}
#endif


	if (working)
	{
	  for (int e = 0; e < surfnr.Size(); e++)
	  {
		  if (surfnr[e] == -1)
		  {
			  continue;
		  }
		  SetThreadPercent((double)e / surfnr.Size() * 100.0);

		  PointIndex pi1 = new PointIndex();
		  PointIndex pi2 = new PointIndex();
		  top.GetEdgeVertices(e+1, pi1, pi2);
		  bool swap = (pi1 > pi2);

		  Point < 3> p1 = mesh[pi1];
		  Point < 3> p2 = mesh[pi2];

		  int order1 = edgeorder[e];
		  int ndof = Math.Max(0, order1 - 1);

		  if (rational && order1 >= 2)
		  {
		  Point < 3> pm = netgen.GlobalMembers.Center(p1, p2);

		  Vec < 3> n1 = @ref.GetNormal(p1, surfnr[e], gi0[e]);
		  Vec < 3> n2 = @ref.GetNormal(p2, surfnr[e], gi1[e]);

		  // p3 = pm + alpha1 n1 + alpha2 n2

		  Mat < 2> mat, inv;
		  Vec < 2> rhs, sol;

		  mat(0,0) = n1 * n1;
		  mat(0,1) = mat(1,0) = n1 * n2;
		  mat(1,1) = n2 * n2;

		  rhs(0) = n1 * (p1 - pm);
		  rhs(1) = n2 * (p2 - pm);


		  Point < 3> p3;

		  if (ngsimd.GlobalMembers.fabs(Det(mat)) > 1e-10)
		  {
			  netgen.GlobalMembers.CalcInverse(mat, inv);
			  sol = inv * rhs;

			  p3 = pm + sol(0) * n1 + sol(1) * n2;
		  }
		  else
		  {
			p3 = pm;
		  }

		  edgecoeffs[edgecoeffsindex[e]] = Vec < 3> (p3);


		  double wold = 1;
		  double w = 1;
		  double dw = 0.1;
		  double dold = 1e99;
		  while (ngsimd.GlobalMembers.fabs(dw) > 1e-12)
		  {
			  Vec < 3> v05 = 0.25 * Vec < 3> (p1) + 0.5 * w * Vec < 3>(p3) + 0.25 * Vec < 3> (p2);
			  v05 /= 1 + (w - 1) * 0.5;
			  Point < 3> p05(v05), pp05(v05);
			  @ref.ProjectToSurface(pp05, surfnr[e], gi0[e]);
			  double d = netgen.GlobalMembers.Dist(pp05, p05);

			  if (d < dold)
			  {
			  dold = d;
			  wold = w;
			  w += dw;
			  }
			  else
			  {
			  dw *= -0.7;
			  w = wold + dw;
			  }
		  }

		  edgeweight[e] = w;
		  continue;
		  }

		  Vector shape = new Vector(ndof);
		  DenseMatrix mat = new DenseMatrix(ndof, ndof);
		  DenseMatrix inv = new DenseMatrix(ndof, ndof);
		  DenseMatrix rhs = new DenseMatrix(ndof, 3);
		  DenseMatrix sol = new DenseMatrix(ndof, 3);

		  rhs = 0.0;
		  mat = 0.0;
		  for (int j = 0; j < xi.Size(); j++)
		  {
		  Point < 3> p;
		  Point < 3> pp;
		  PointGeomInfo ppgi = new PointGeomInfo();

		  if (swap)
		  {
			  p = p1 + xi[j] * (p2 - p1);
			  @ref.PointBetween(p1, p2, xi[j], surfnr[e], gi0[e], gi1[e], pp, ppgi);
		  }
		  else
		  {
			  p = p2 + xi[j] * (p1 - p2);
			  @ref.PointBetween(p2, p1, xi[j], surfnr[e], gi1[e], gi0[e], pp, ppgi);
		  }

		  Vec < 3> dist = pp - p;

		  netgen.GlobalMembers.CalcEdgeShape(order1, 2 * xi[j] - 1, shape(0));

		  for (int k = 0; k < ndof; k++)
		  {
			for (int l = 0; l < ndof; l++)
			{
			  mat(k,l) += weight[j] * shape(k) * shape(l);
			}
		  }

		  for (int k = 0; k < ndof; k++)
		  {
			for (int l = 0; l < 3; l++)
			{
			  rhs(k,l) += weight[j] * shape(k) * dist(l);
			}
		  }
		  }

		  netgen.GlobalMembers.CalcInverse(mat, inv);
		  Mult(inv, rhs, sol);

		  int first = edgecoeffsindex[e];
		  for (int j = 0; j < ndof; j++)
		  {
		for (int k = 0; k < 3; k++)
		{
		  edgecoeffs[first + j](k) = sol(j,k);
		}
		  }
	  }
	}
	}


	Array<int> use_edge = new Array<int>(nedges);
	Array<int> edge_surfnr1 = new Array<int>(nedges);
	Array<int> edge_surfnr2 = new Array<int>(nedges);
	Array<int> swap_edge = new Array<int>(nedges);
	Array<EdgePointGeomInfo> edge_gi0 = new Array<EdgePointGeomInfo>(nedges);
	Array<EdgePointGeomInfo> edge_gi1 = new Array<EdgePointGeomInfo>(nedges);
	use_edge = 0;

	if (working)
	{
	  for (SegmentIndex i = 0; i < mesh.GetNSeg(); i++)
	  {
	  Segment seg = mesh[i];
	  int edgenr = top.GetEdge(i);
	  use_edge[edgenr] = 1;
	  edge_surfnr1[edgenr] = seg.surfnr1;
	  edge_surfnr2[edgenr] = seg.surfnr2;
	  edge_gi0[edgenr] = seg.epgeominfo[0];
	  edge_gi1[edgenr] = seg.epgeominfo[1];
	  swap_edge[edgenr] = (int)(seg[0] > seg[1]);
	  }
	}

#if PARALLEL
	if (ntasks > 1)
	{
	// distribute it ...
	TABLE<double> senddata = new TABLE<double>(ntasks);
	TABLE<double> recvdata = new TABLE<double>(ntasks);
	if (working)
	{
	  for (int e = 0; e < nedges; e++)
	  {
		  partop.GetDistantEdgeNums(e+1, procs);
		  for (int j = 0; j < procs.Size(); j++)
		  {
		  senddata.Add(procs[j], use_edge[e]);
		  if (use_edge[e])
		  {
			  senddata.Add(procs[j], edge_surfnr1[e]);
			  senddata.Add(procs[j], edge_surfnr2[e]);
			  senddata.Add(procs[j], edge_gi0[e].edgenr);
			  senddata.Add(procs[j], edge_gi0[e].body);
			  senddata.Add(procs[j], edge_gi0[e].dist);
			  senddata.Add(procs[j], edge_gi0[e].u);
			  senddata.Add(procs[j], edge_gi0[e].v);
			  senddata.Add(procs[j], edge_gi1[e].edgenr);
			  senddata.Add(procs[j], edge_gi1[e].body);
			  senddata.Add(procs[j], edge_gi1[e].dist);
			  senddata.Add(procs[j], edge_gi1[e].u);
			  senddata.Add(procs[j], edge_gi1[e].v);
			  senddata.Add(procs[j], swap_edge[e]);
		  }
		  }
	  }
	}
	netgen.GlobalMembers.MyMPI_ExchangeTable(senddata, recvdata, MPI_TAG_CURVE, curve_comm);
	Array<int> cnt = new Array<int>(ntasks);
	cnt = 0;
	if (working)
	{
	  for (int e = 0; e < edge_surfnr1.Size(); e++)
	  {
		  partop.GetDistantEdgeNums(e+1, procs);
		  for (int j = 0; j < procs.Size(); j++)
		  {
		  int get_edge = recvdata[procs[j]][cnt[procs[j]]++];
		  if (get_edge != 0)
		  {
			  use_edge[e] = 1;
			  edge_surfnr1[e] = (int)recvdata[procs[j]][cnt[procs[j]]++];
			  edge_surfnr2[e] = (int)recvdata[procs[j]][cnt[procs[j]]++];
			  edge_gi0[e].edgenr = (int)recvdata[procs[j]][cnt[procs[j]]++];
			  edge_gi0[e].body = (int)recvdata[procs[j]][cnt[procs[j]]++];
			  edge_gi0[e].dist = recvdata[procs[j]][cnt[procs[j]]++];
			  edge_gi0[e].u = recvdata[procs[j]][cnt[procs[j]]++];
			  edge_gi0[e].v = recvdata[procs[j]][cnt[procs[j]]++];
			  edge_gi1[e].edgenr = (int)recvdata[procs[j]][cnt[procs[j]]++];
			  edge_gi1[e].body = (int)recvdata[procs[j]][cnt[procs[j]]++];
			  edge_gi1[e].dist = recvdata[procs[j]][cnt[procs[j]]++];
			  edge_gi1[e].u = recvdata[procs[j]][cnt[procs[j]]++];
			  edge_gi1[e].v = recvdata[procs[j]][cnt[procs[j]]++];
			  swap_edge[e] = recvdata[procs[j]][cnt[procs[j]]++];
		  }
		  }
	  }
	}

	}
#endif

	if (working)
	{
	  for (int edgenr = 0; edgenr < use_edge.Size(); edgenr++)
	  {
	  int segnr = edgenr;
	  if (!use_edge[edgenr])
	  {
		  continue;
	  }

	  SetThreadPercent((double)edgenr / edge_surfnr1.Size() * 100.0);

	  PointIndex pi1 = new PointIndex();
	  PointIndex pi2 = new PointIndex();
	  top.GetEdgeVertices(edgenr + 1, pi1, pi2);

	  bool swap = swap_edge[edgenr]; // (pi1 > pi2);
	  if (swap)
	  {
		  netgen.GlobalMembers.Swap(ref pi1, ref pi2);
	  }

	  Point < 3> p1 = mesh[pi1];
	  Point < 3> p2 = mesh[pi2];

	  int order1 = edgeorder[segnr];
	  int ndof = Math.Max(0, order1 - 1);

	  if (rational)
	  {
		  Vec < 3> tau1 = @ref.GetTangent(p1, edge_surfnr2[edgenr], edge_surfnr1[edgenr], edge_gi0[edgenr]);
		  Vec < 3> tau2 = @ref.GetTangent(p2, edge_surfnr2[edgenr], edge_surfnr1[edgenr], edge_gi1[edgenr]);
		  // p1 + alpha1 tau1 = p2 + alpha2 tau2;

		  Mat < 3,2> mat;
		  Mat < 2,3> inv;
		  Vec < 3> rhs;
		  Vec < 2> sol;
		  for (int j = 0; j < 3; j++)
		  {
		  mat(j,0) = tau1(j);
		  mat(j,1) = -tau2(j);
		  rhs(j) = p2(j) - p1(j);
		  }
		  netgen.GlobalMembers.CalcInverse(mat, inv);
		  sol = inv * rhs;

		  Point < 3> p3 = p1 + sol(0) * tau1;
		  edgecoeffs[edgecoeffsindex[segnr]] = Vec < 3> (p3);

		  double wold = 1;
		  double w = 1;
		  double dw = 0.1;
		  double dold = 1e99;
		  while (ngsimd.GlobalMembers.fabs(dw) > 1e-12)
		  {
		  Vec < 3> v05 = 0.25 * Vec < 3> (p1) + 0.5 * w * Vec < 3>(p3) + 0.25 * Vec < 3> (p2);
		  v05 /= 1 + (w - 1) * 0.5;
		  Point < 3> p05(v05), pp05(v05);
		  @ref.ProjectToEdge(pp05, edge_surfnr1[edgenr], edge_surfnr2[edgenr], edge_gi0[edgenr]);
		  double d = netgen.GlobalMembers.Dist(pp05, p05);

		  if (d < dold)
		  {
			  dold = d;
			  wold = w;
			  w += dw;
		  }
		  else
		  {
			  dw *= -0.7;
			  w = wold + dw;
		  }
		  // *testout << "w = " << w << ", dw = " << dw << endl;
		  }

		  // cout << "wopt = " << w << ", dopt = " << dold << endl;
		  edgeweight[segnr] = w;

		  //             cout << "p1 = " << p1 << ", tau1 = " << tau1 << ", alpha1 = " << sol(0) << endl;
		  //             cout << "p2 = " << p2 << ", tau2 = " << tau2 << ", alpha2 = " << -sol(1) << endl;
		  //             cout << "p+alpha tau = " << p1 + sol(0) * tau1
		  //                  << " =?= " << p2 +sol(1) * tau2 << endl;

	  }

	  else

	  {
		  Vector shape = new Vector(ndof);
		  DenseMatrix mat = new DenseMatrix(ndof, ndof);
		  DenseMatrix inv = new DenseMatrix(ndof, ndof);
		  DenseMatrix rhs = new DenseMatrix(ndof, 3);
		  DenseMatrix sol = new DenseMatrix(ndof, 3);

		  rhs = 0.0;
		  mat = 0.0;
		  for (int j = 0; j < xi.Size(); j++)
		  {
		  Point < 3> p, pp;
		  EdgePointGeomInfo ppgi = new EdgePointGeomInfo();

		  if (swap)
		  {
			  p = p1 + xi[j] * (p2 - p1);
			  @ref.PointBetween(p1, p2, xi[j], edge_surfnr2[edgenr], edge_surfnr1[edgenr], edge_gi0[edgenr], edge_gi1[edgenr], pp, ppgi);
		  }
		  else
		  {
			  p = p2 + xi[j] * (p1 - p2);
			  @ref.PointBetween(p2, p1, xi[j], edge_surfnr2[edgenr], edge_surfnr1[edgenr], edge_gi1[edgenr], edge_gi0[edgenr], pp, ppgi);
		  }

		  Vec < 3> dist = pp - p;

		  netgen.GlobalMembers.CalcEdgeShape(order1, 2 * xi[j] - 1, shape(0));

		  for (int k = 0; k < ndof; k++)
		  {
			for (int l = 0; l < ndof; l++)
			{
			  mat(k,l) += weight[j] * shape(k) * shape(l);
			}
		  }

		  for (int k = 0; k < ndof; k++)
		  {
			for (int l = 0; l < 3; l++)
			{
			  rhs(k,l) += weight[j] * shape(k) * dist(l);
			}
		  }
		  }


		  netgen.GlobalMembers.CalcInverse(mat, inv);
		  Mult(inv, rhs, sol);

		  int first = edgecoeffsindex[segnr];
		  for (int j = 0; j < ndof; j++)
		  {
		for (int k = 0; k < 3; k++)
		{
		  edgecoeffs[first + j](k) = sol(j,k);
		}
		  }
	  }
	  }
	}



	PrintMessage(3, "Curving faces");

	Array<int> surfnr = new Array<int>(nfaces);
	surfnr = -1;

	if (working)
	{
	  for (SurfaceElementIndex i = 0; i < mesh.GetNSE(); i++)
	  {
	surfnr[top.GetFace(i)] = mesh.GetFaceDescriptor(mesh[i].GetIndex()).SurfNr();
	  }
	}

#if PARALLEL
	TABLE<int> send_surfnr = new TABLE<int>(ntasks);
	TABLE<int> recv_surfnr = new TABLE<int>(ntasks);

	if (ntasks > 1 && working)
	{
	for (int f = 0; f < nfaces; f++)
	{
		partop.GetDistantFaceNums(f + 1, procs);
		for (int j = 0; j < procs.Size(); j++)
		{
		  send_surfnr.Add(procs[j], surfnr[f]);
		}
	}
	}

	if (ntasks > 1)
	{
	  netgen.GlobalMembers.MyMPI_ExchangeTable(send_surfnr, recv_surfnr, MPI_TAG_CURVE, curve_comm);
	}

	if (ntasks > 1 && working)
	{
	Array<int> cnt = new Array<int>(ntasks);
	cnt = 0;
	for (int f = 0; f < nfaces; f++)
	{
		partop.GetDistantFaceNums(f + 1, procs);
		for (int j = 0; j < procs.Size(); j++)
		{
		  surfnr[f] = Math.Max(surfnr[f], recv_surfnr[procs[j]][cnt[procs[j]]++]);
		}
	}
	}
#endif

	if (mesh.GetDimension() == 3 && working)
	{
	for (int f = 0; f < nfaces; f++)
	{
		int facenr = f;
		if (surfnr[f] == -1)
		{
			continue;
		}
		// if (el.GetType() == TRIG && order >= 3)
		if (top.GetFaceType(facenr + 1) == TRIG && order >= 3)
		{
		ArrayMem<int, 3> verts = new ArrayMem<int, 3>(3);
		top.GetFaceVertices(facenr + 1, verts);

		int[] fnums = {0, 1, 2};
		/*
		if (el[fnums[0]] > el[fnums[1]]) swap (fnums[0], fnums[1]);
		if (el[fnums[1]] > el[fnums[2]]) swap (fnums[1], fnums[2]);
		if (el[fnums[0]] > el[fnums[1]]) swap (fnums[0], fnums[1]);
		*/
		if (verts[fnums[0]] > verts[fnums[1]])
		{
			swap(fnums[0], fnums[1]);
		}
		if (verts[fnums[1]] > verts[fnums[2]])
		{
			swap(fnums[1], fnums[2]);
		}
		if (verts[fnums[0]] > verts[fnums[1]])
		{
			swap(fnums[0], fnums[1]);
		}

		int order1 = faceorder[facenr];
		int ndof = Math.Max(0, (order1 - 1) * (order1 - 2) / 2);

		Vector shape = new Vector(ndof);
		Vector dmat = new Vector(ndof);
		MatrixFixWidth < 3> rhs(ndof), sol(ndof);

		rhs = 0.0;
		dmat = 0.0;

		int np = netgen.GlobalMembers.sqr(xi.Size());
		Array<Point < 2>> xia = new Array<Point < 2>>(np);
		Array<Point < 3>> xa = new Array<Point < 3>>(np);

		for (int jx = 0, jj = 0; jx < xi.Size(); jx++)
		{
		  for (int jy = 0; jy < xi.Size(); jy++, jj++)
		  {
			xia[jj] = Point < 2> ((1 - xi[jy]) * xi[jx], xi[jy]);
		  }
		}

		// CalcMultiPointSurfaceTransformation (&xia, i, &xa, NULL);

		Array<int> edgenrs = new Array<int>();
		top.GetFaceEdges(facenr + 1, edgenrs);
		for (int k = 0; k < edgenrs.Size(); k++)
		{
			edgenrs[k]--;
		}

		for (int jj = 0; jj < np; jj++)
		{
			Point < 3> pp(0,0,0);
			double[] lami = {xia[jj](0), xia[jj](1), 1 - xia[jj](0) - xia[jj](1)};

			for (int k = 0; k < verts.Size(); k++)
			{
			  pp += lami[k] * Vec < 3> (new mesh.Point(verts[k]));
			}

			// const ELEMENT_EDGE * edges = MeshTopology::GetEdges0 (TRIG);
			for (int k = 0; k < edgenrs.Size(); k++)
			{
			int eorder = edgeorder[edgenrs[k]];
			if (eorder < 2)
			{
				continue;
			}

			int first = edgecoeffsindex[edgenrs[k]];
			Vector eshape = new Vector(eorder - 1);
			int vi1;
			int vi2;
			top.GetEdgeVertices(edgenrs[k] + 1, vi1, vi2);
			if (vi1 > vi2)
			{
				swap(vi1, vi2);
			}
			int v1 = -1;
			int v2 = -1;
			for (int j = 0; j < 3; j++)
			{
				if (verts[j] == vi1)
				{
					v1 = j;
				}
				if (verts[j] == vi2)
				{
					v2 = j;
				}
			}

			netgen.GlobalMembers.CalcScaledEdgeShape(eorder, lami[v1] - lami[v2], lami[v1] + lami[v2], eshape(0));
			for (int n = 0; n < eshape.Size(); n++)
			{
			  pp += eshape(n) * edgecoeffs[first + n];
			}
			}
			xa[jj] = pp;
		}

		for (int jx = 0, jj = 0; jx < xi.Size(); jx++)
		{
		  for (int jy = 0; jy < xi.Size(); jy++, jj++)
		  {
			  double y = xi[jy];
			  double x = (1 - y) * xi[jx];
			  double[] lami = {x, y, 1 - x - y};
			  double wi = weight[jx] * weight[jy] * (1 - y);

			  Point < 3> pp = xa[jj];
			  // ref -> ProjectToSurface (pp, mesh.GetFaceDescriptor(el.GetIndex()).SurfNr());
			  @ref.ProjectToSurface(pp, surfnr[facenr]);
			  Vec < 3> dist = pp - xa[jj];

			  netgen.GlobalMembers.CalcTrigShape(order1, lami[fnums[1]] - lami[fnums[0]], 1 - lami[fnums[1]] - lami[fnums[0]], shape(0));

			  for (int k = 0; k < ndof; k++)
			  {
			dmat(k) += wi * shape(k) * shape(k);
			  }

			  dist *= wi;
			  for (int k = 0; k < ndof; k++)
			  {
			for (int l = 0; l < 3; l++)
			{
			  rhs(k,l) += shape(k) * dist(l);
			}
			  }
		  }
		}

		for (int i = 0; i < ndof; i++)
		{
		  for (int j = 0; j < 3; j++)
		  {
			sol(i,j) = rhs(i,j) / dmat(i); // Orthogonal basis !
		  }
		}

		int first = facecoeffsindex[facenr];
		for (int j = 0; j < ndof; j++)
		{
		  for (int k = 0; k < 3; k++)
		  {
			facecoeffs[first + j](k) = sol(j,k);
		  }
		}
		}
	}
	}


	// compress edge and face tables
	int newbase = 0;
	for (int i = 0; i < edgeorder.Size(); i++)
	{
	bool curved = false;
	int oldbase = edgecoeffsindex[i];
	int nd = edgecoeffsindex[i + 1] - edgecoeffsindex[i];

	for (int j = 0; j < nd; j++)
	{
	  if (edgecoeffs[oldbase + j].Length() > 1e-12)
	  {
		curved = true;
	  }
	}
	if (rational)
	{
		curved = true;
	}

	if (curved && newbase != oldbase)
	{
	  for (int j = 0; j < nd; j++)
	  {
		edgecoeffs[newbase + j] = edgecoeffs[oldbase + j];
	  }
	}

	edgecoeffsindex[i] = newbase;
	if (!curved)
	{
		edgeorder[i] = 1;
	}
	if (curved)
	{
		newbase += nd;
	}
	}
	edgecoeffsindex.Last() = newbase;


	newbase = 0;
	for (int i = 0; i < faceorder.Size(); i++)
	{
	bool curved = false;
	int oldbase = facecoeffsindex[i];
	int nd = facecoeffsindex[i + 1] - facecoeffsindex[i];

	for (int j = 0; j < nd; j++)
	{
	  if (facecoeffs[oldbase + j].Length() > 1e-12)
	  {
		curved = true;
	  }
	}

	if (curved && newbase != oldbase)
	{
	  for (int j = 0; j < nd; j++)
	  {
		facecoeffs[newbase + j] = facecoeffs[oldbase + j];
	  }
	}

	facecoeffsindex[i] = newbase;
	if (!curved)
	{
		faceorder[i] = 1;
	}
	if (curved)
	{
		newbase += nd;
	}
	}
	facecoeffsindex.Last() = newbase;

	if (working)
	{
	  ishighorder = (order > 1);
	}
	// (*testout) << "edgecoeffs = " << endl << edgecoeffs << endl;
	// (*testout) << "facecoeffs = " << endl << facecoeffs << endl;


#if PARALLEL
	curve_comm.Barrier();
	// MPI_Comm_free (&curve_comm);
#endif
  }

  public int GetOrder()
  {
	  return order;
  }

  public virtual void DoArchive(Archive ar)
  {
	ar & edgeorder & faceorder & edgecoeffsindex & facecoeffsindex & edgecoeffs & facecoeffs & edgeweight & order & rational & ishighorder;
  }


  // ***********************  Transform edges *****************************


//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool IsSegmentCurved(SegmentIndex elnr) const
  public bool IsSegmentCurved(SegmentIndex elnr)
  {
	if (mesh.coarsemesh)
	{
	HPRefElement hpref_el = (*mesh.hpelements)[mesh[elnr].hp_elnr];

	return mesh.coarsemesh.GetCurvedElements().IsSegmentCurved(hpref_el.coarse_elnr);
	}

	SegmentInfo info = new SegmentInfo();
	info.elnr = elnr;
	info.order = order;
	info.ndof = info.nv = 2;
	if (info.order > 1)
	{
	MeshTopology top = mesh.GetTopology();
	info.edgenr = top.GetSegmentEdge(elnr + 1) - 1;
	info.ndof += edgeorder[info.edgenr] - 1;
	}

	return (info.ndof > info.nv);
  }


  // ********************** Transform surface elements *******************


//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool IsSurfaceElementCurved(SurfaceElementIndex elnr) const
  public bool IsSurfaceElementCurved(SurfaceElementIndex elnr)
  {
	if (mesh[elnr].GetType() != TRIG)
	{
		return true;
	}
	if (!IsHighOrder())
	{
		return false;
	}

	if (mesh.coarsemesh)
	{
	HPRefElement hpref_el = (*mesh.hpelements)[mesh[elnr].hp_elnr];

	return mesh.coarsemesh.GetCurvedElements().IsSurfaceElementCurved(hpref_el.coarse_elnr);
	}

	Element2d el = mesh[elnr];
	ELEMENT_TYPE type = el.GetType();

	SurfaceElementInfo info = new SurfaceElementInfo();
	info.elnr = elnr;
	info.order = order;

	switch (type)
	{
	  case TRIG :
		  info.nv = 3;
		  break;
	  case QUAD :
		  info.nv = 4;
		  break;
	  case TRIG6:
		  return true;
	  default:
	cerr << "undef element in CalcSurfaceTrafo" << "\n";
	break;
	}
	info.ndof = info.nv;

	// info.ndof = info.nv = ( (type == TRIG) || (type == TRIG6) ) ? 3 : 4;
	if (info.order > 1)
	{
	MeshTopology top = mesh.GetTopology();

	top.GetSurfaceElementEdges(elnr + 1, info.edgenrs);
	for (int i = 0; i < info.edgenrs.Size(); i++)
	{
	  info.edgenrs[i]--;
	}
	info.facenr = top.GetSurfaceElementFace(elnr + 1) - 1;

	for (int i = 0; i < info.edgenrs.Size(); i++)
	{
	  info.ndof += edgecoeffsindex[info.edgenrs[i] + 1] - edgecoeffsindex[info.edgenrs[i]];
	}
	info.ndof += facecoeffsindex[info.facenr + 1] - facecoeffsindex[info.facenr];
	}

	return (info.ndof > info.nv);
  }


  // ********************** Transform volume elements *******************


//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool IsElementCurved(ElementIndex elnr) const
  public bool IsElementCurved(ElementIndex elnr)
  {
	if (mesh[elnr].GetType() != TET)
	{
		return true;
	}

	if (mesh.coarsemesh)
	{
	HPRefElement hpref_el = (*mesh.hpelements)[mesh[elnr].hp_elnr];

	return mesh.coarsemesh.GetCurvedElements().IsElementCurved(hpref_el.coarse_elnr);
	}

	Element el = mesh[elnr];
	ELEMENT_TYPE type = el.GetType();

	int nfaces = MeshTopology.GetNFaces(type);
	if (nfaces > 4)
	{ // not a tet
	ELEMENT_FACE faces = MeshTopology.GetFaces0(type);
	for (int j = 0; j < nfaces; j++)
	{
		if (faces[j][3] != -1)
		{ // a quad face
		Point < 3> pts[4];
		for (int k = 0; k < 4; k++)
		{
		  pts[k] = new mesh.Point(el[faces[j][k]]);
		}
		Vec < 3> twist = (pts[1] - pts[0]) - (pts[2] - pts[3]);
		if (twist.Length() > 1e-8 * (pts[1] - pts[0]).Length())
		{
		  return true;
		}
		}
	}
	}



	ElementInfo info = new ElementInfo();
	info.elnr = elnr;
	info.order = order;
	info.ndof = info.nv = MeshTopology.GetNPoints(type);
	if (info.order > 1)
	{
	MeshTopology top = mesh.GetTopology();

	info.nedges = top.GetElementEdges(elnr + 1, info.edgenrs, 0);
	for (int i = 0; i < info.nedges; i++)
	{
	  info.edgenrs[i]--;
	}

	info.nfaces = top.GetElementFaces(elnr + 1, info.facenrs, 0);
	for (int i = 0; i < info.nfaces; i++)
	{
	  info.facenrs[i]--;
	}

	for (int i = 0; i < info.nedges; i++)
	{
	  info.ndof += edgecoeffsindex[info.edgenrs[i] + 1] - edgecoeffsindex[info.edgenrs[i]];
	}
	for (int i = 0; i < info.nfaces; i++)
	{
	  info.ndof += facecoeffsindex[info.facenrs[i] + 1] - facecoeffsindex[info.facenrs[i]];
	}
	}

	return (info.ndof > info.nv);
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool IsElementHighOrder(ElementIndex elnr) const
  public bool IsElementHighOrder(ElementIndex elnr)
  {
	if (mesh.coarsemesh)
	{
	HPRefElement hpref_el = (*mesh.hpelements)[mesh[elnr].hp_elnr];

	return mesh.coarsemesh.GetCurvedElements().IsElementHighOrder(hpref_el.coarse_elnr);
	}

	Element el = mesh[elnr];
	ELEMENT_TYPE type = el.GetType();

	ElementInfo info = new ElementInfo();
	info.elnr = elnr;
	info.order = order;
	info.ndof = info.nv = MeshTopology.GetNPoints(type);
	if (info.order > 1)
	{
	MeshTopology top = mesh.GetTopology();

	info.nedges = top.GetElementEdges(elnr + 1, info.edgenrs, 0);
	for (int i = 0; i < info.nedges; i++)
	{
		info.edgenrs[i]--;
	}

	info.nfaces = top.GetElementFaces(elnr + 1, info.facenrs, 0);
	for (int i = 0; i < info.nfaces; i++)
	{
		info.facenrs[i]--;
	}

	for (int i = 0; i < info.nedges; i++)
	{
		  if (edgecoeffsindex[info.edgenrs[i] + 1] > edgecoeffsindex[info.edgenrs[i]])
		  {
			  return true;
		  }
	}
	for (int i = 0; i < info.nfaces; i++)
	{
		  if (facecoeffsindex[info.facenrs[i] + 1] > facecoeffsindex[info.facenrs[i]])
		  {
			  return true;
		  }
	}
	}
	return false;
  }


  public void CalcSegmentTransformation(double xi, SegmentIndex segnr, Point < 3> x)
  {
	  CalcSegmentTransformation<double> (xi, segnr, x, null);
  }

  public void CalcSegmentTransformation(double xi, SegmentIndex segnr, Vec < 3> dxdxi)
  {
	  CalcSegmentTransformation<double> (xi, segnr, null, dxdxi);
  }

  public void CalcSegmentTransformation(double xi, SegmentIndex segnr, Point < 3> x, Vec < 3> dxdxi)
  {
	  CalcSegmentTransformation<double> (xi, segnr, x, dxdxi, null);
  }

  public void CalcSegmentTransformation(double xi, SegmentIndex segnr, Point < 3> x, Vec < 3> dxdxi, ref bool curved)
  {
	  CalcSegmentTransformation(xi, segnr, x, dxdxi, curved);
  }



  public void CalcSurfaceTransformation(Point < 2> xi, SurfaceElementIndex elnr, Point < 3> x)
  {
	  CalcSurfaceTransformation(xi, elnr, x, null);
  }

  public void CalcSurfaceTransformation(Point < 2> xi, SurfaceElementIndex elnr, Mat<3, ref 2> dxdxi)
  {
	  CalcSurfaceTransformation(xi, elnr, null, dxdxi);
  }

  public void CalcSurfaceTransformation(Point < 2> xi, SurfaceElementIndex elnr, Point < 3> x, Mat<3, ref 2> dxdxi)
  {
	  CalcSurfaceTransformation(xi, elnr, x, dxdxi, null);
  }

  public void CalcSurfaceTransformation(Point < 2> xi, SurfaceElementIndex elnr, Point < 3> x, Mat<3, ref 2> dxdxi, ref bool curved)
  {
	  CalcSurfaceTransformation(xi, elnr, x, dxdxi, curved);
  }





  public void CalcElementTransformation(Point < 3> xi, ElementIndex elnr, Point < 3> x)
  {
	  CalcElementTransformation(xi, elnr, x, null);
  }

  public void CalcElementTransformation(Point < 3> xi, ElementIndex elnr, Mat<3, ref 3> dxdxi)
  {
	  CalcElementTransformation(xi, elnr, null, dxdxi);
  }

  public void CalcElementTransformation(Point < 3> xi, ElementIndex elnr, Point < 3> x, Mat<3, ref 3> dxdxi)
  {
	  CalcElementTransformation(xi, elnr, x, dxdxi);
  }

  public void CalcElementTransformation(Point < 3> xi, ElementIndex elnr, Point < 3> x, Mat<3, ref 3> dxdxi, object buffer, bool valid)
  {
	  CalcElementTransformation(xi, elnr, x, dxdxi, buffer, valid);
  }

  // void CalcElementTransformation (const Point<3> & xi, ElementIndex elnr,
  // 				  Point<3> & x, Mat<3,3> & dxdxi) // , bool & curved)
  //   { CalcElementTransformation (xi, elnr, &x, &dxdxi /* , &curved * ); }


  /*
  void CalcMultiPointSegmentTransformation (Array<double> * xi, SegmentIndex segnr,
						Array<Point<3> > * x,
						Array<Vec<3> > * dxdxi);
  */

//C++ TO C# CONVERTER TODO TASK: C++ template specifiers with non-type parameters cannot be converted to C#:
//ORIGINAL LINE: template <int DIM_SPACE, typename T>

  /*
  void CurvedElements :: 
  CalcMultiPointSegmentTransformation (Array<double> * xi, SegmentIndex segnr,
					   Array<Point<3> > * x,
					   Array<Vec<3> > * dxdxi)
  {
    ;
  }
  */

  public void CalcMultiPointSegmentTransformation<int DIM_SPACE, T>(SegmentIndex elnr, int n, T[] xi, uint sxi, T[] x, uint sx, T[] dxdxi, uint sdxdxi)
  {
	for (int ip = 0; ip < n; ip++)
	{
	Point<3,T> xg = new Point<3,T>();
	Vec<3,T> dx = new Vec<3,T>();

	// mesh->GetCurvedElements().
	CalcSegmentTransformation<T> (xi[ip * sxi], elnr, xg, dx);

	if (x != null)
	{
	  for (int i = 0; i < DIM_SPACE; i++)
	  {
		x[ip * sx + i] = xg(i);
	  }
	}

	if (dxdxi != null)
	{
	  for (int i = 0; i < DIM_SPACE; i++)
	  {
		dxdxi[ip * sdxdxi + i] = dx(i);
	  }
	}
	}
  }

  public void CalcMultiPointSurfaceTransformation(Array< Point < 2>> xi, SurfaceElementIndex elnr, Array< Point < 3>> x, Array< Mat < 3,2>> dxdxi)
  {
//C++ TO C# CONVERTER TODO TASK: C# does not have an equivalent to pointers to value types:
//ORIGINAL LINE: double * px = (x) ? &(*x)[0](0) : null;
	double px = (x) != null ?  x[0](0) : null;
//C++ TO C# CONVERTER TODO TASK: C# does not have an equivalent to pointers to value types:
//ORIGINAL LINE: double * pdxdxi = (dxdxi) ? &(*dxdxi)[0](0) : null;
	double pdxdxi = (dxdxi) != null ?  dxdxi[0](0) : null;

	CalcMultiPointSurfaceTransformation < 3> (elnr, xi.Size(), xi[0](0), 2, px, 3, pdxdxi, 6);
  }

//C++ TO C# CONVERTER TODO TASK: C++ template specifiers with non-type parameters cannot be converted to C#:
//ORIGINAL LINE: template <int DIM_SPACE, typename T>
  public void CalcMultiPointSurfaceTransformation<int DIM_SPACE, T>(SurfaceElementIndex elnr, int npts, T[] xi, uint sxi, T[] x, uint sx, T[] dxdxi, uint sdxdxi)
  {
	if (mesh.coarsemesh)
	{
	HPRefElement hpref_el = (*mesh.hpelements)[mesh[elnr].hp_elnr];

	// xi umrechnen
	T[] lami = Arrays.InitializeWithDefaultInstances<T>(4);
	TFlatVector<T> vlami = new TFlatVector<T>(4, lami);

	ArrayMem<Point<2,T>, 50> coarse_xi = new ArrayMem<Point<2,T>, 50>((uint)npts);

	for (int pi = 0; pi < npts; pi++)
	{
		vlami = 0;
		Point<2,T> hxi = new Point<2,T>(xi[pi * sxi], xi[pi * sxi + 1]);
		mesh[elnr].GetShapeNew(hxi, vlami);

		Point<2,T> cxi = new Point<2,T>(0, 0);
		for (int i = 0; i < hpref_el.np; i++)
		{
		  for (int j = 0; j < 2; j++)
		  {
		cxi(j) += hpref_el.param[i][j] * lami[i];
		  }
		}

		coarse_xi[pi] = cxi;
	}

	mesh.coarsemesh.GetCurvedElements().CalcMultiPointSurfaceTransformation<DIM_SPACE,T> (hpref_el.coarse_elnr, npts, coarse_xi[0](0), coarse_xi[1](0) - coarse_xi[0](0), x, sx, dxdxi, sdxdxi);

	// Mat<3,2> dxdxic;
	if (dxdxi != null)
	{
			T[] mem_dlami = Arrays.InitializeWithDefaultInstances<T>(8); // avoid alignment problems if T is SIMD
		MatrixFixWidth<2,T> dlami = new MatrixFixWidth<2,T>(4, mem_dlami);
		dlami = T(0.0);

		for (int pi = 0; pi < npts; pi++)
		{
		Point<2,T> hxi = new Point<2,T>(xi[pi * sxi], xi[pi * sxi + 1]);
		mesh[elnr].GetDShapeNew(hxi, dlami);

		Mat<2,2,T> trans = new Mat<2,2,T>();
		trans = 0;
		for (int k = 0; k < 2; k++)
		{
		  for (int l = 0; l < 2; l++)
		  {
			for (int i = 0; i < hpref_el.np; i++)
			{
			  trans(l,k) += hpref_el.param[i][l] * dlami(i, k);
			}
		  }
		}

		Mat<DIM_SPACE,2,T> hdxdxic = new Mat<DIM_SPACE,2,T>();
		Mat<DIM_SPACE,2,T> hdxdxi = new Mat<DIM_SPACE,2,T>();
		for (int k = 0; k < 2 * DIM_SPACE; k++)
		{
		  hdxdxic(k) = dxdxi[pi * sdxdxi + k];
		}

		hdxdxi = hdxdxic * trans;

		for (int k = 0; k < 2 * DIM_SPACE; k++)
		{
		  dxdxi[pi * sdxdxi + k] = hdxdxi(k);
		}

		// dxdxic = (*dxdxi)[pi];
		// (*dxdxi)[pi] = dxdxic * trans;
		}
	}

	return;
	}


	Element2d el = mesh[elnr];
	ELEMENT_TYPE type = el.GetType();

	SurfaceElementInfo info = new SurfaceElementInfo();
	info.elnr = elnr;
	info.order = order;
	switch (type)
	{
	  case TRIG :
		  info.nv = 3;
		  break;
	  case QUAD :
		  info.nv = 4;
		  break;
	  case TRIG6:
		  info.nv = 6;
		  break;
	  case QUAD8 :
		  info.nv = 8;
		  break;
	  default:
	cerr << "undef element in CalcMultPointSurfaceTrafo" << "\n";
	break;
	}
	info.ndof = info.nv;

	// if (info.order > 1)
	//   {
	//     const MeshTopology & top = mesh.GetTopology();

	//     top.GetSurfaceElementEdges (elnr+1, info.edgenrs);
	//     for (int i = 0; i < info.edgenrs.Size(); i++)
	//       info.edgenrs[i]--;
	//     info.facenr = top.GetSurfaceElementFace (elnr+1)-1;

	//     for (int i = 0; i < info.edgenrs.Size(); i++)
	//       info.ndof += edgecoeffsindex[info.edgenrs[i]+1] - edgecoeffsindex[info.edgenrs[i]];
	//     info.ndof += facecoeffsindex[info.facenr+1] - facecoeffsindex[info.facenr];
	//   }

// Michael Woopen: THESE FOLLOWING LINES ARE COPIED FROM CurvedElements::CalcSurfaceTransformation

	if (info.order > 1)
	{
	MeshTopology top = mesh.GetTopology();

	top.GetSurfaceElementEdges(elnr + 1, info.edgenrs);
	for (int i = 0; i < info.edgenrs.Size(); i++)
	{
	  info.edgenrs[i]--;
	}
	info.facenr = top.GetSurfaceElementFace(elnr + 1) - 1;


	bool firsttry = true;
	bool problem = false;

	while (firsttry || problem)
	{
		problem = false;

		for (int i = 0; !problem && i < info.edgenrs.Size(); i++)
		{
		if (info.edgenrs[i] + 1 >= edgecoeffsindex.Size())
		{
		  problem = true;
		}
		else
		{
		  info.ndof += edgecoeffsindex[info.edgenrs[i] + 1] - edgecoeffsindex[info.edgenrs[i]];
		}
		}
		if (info.facenr + 1 >= facecoeffsindex.Size())
		{
		  problem = true;
		}
		else
		{
		  info.ndof += facecoeffsindex[info.facenr + 1] - facecoeffsindex[info.facenr];
		}

		if (problem && !firsttry)
		{
		  throw new Exception("something wrong with curved elements");
		}

		if (problem)
		{
		  BuildCurvedElements(null,order,rational);
		}

		firsttry = false;
	}
	}



	bool ok = true;
	for (int i = 0; i < npts; i++)
	{
		Point<2,T> _xi = new Point<2,T>(xi[i * sxi], xi[i * sxi + 1]);
		Point<DIM_SPACE,T> _x = new Point<DIM_SPACE,T>();
		Mat<DIM_SPACE,2,T> _dxdxi = new Mat<DIM_SPACE,2,T>();
		if (!EvaluateMapping(info, _xi, _x, _dxdxi))
		{
			  ok = false;
			  break;
		}
		// *testout << "x = " << _x << ", dxdxi = " << _dxdxi << endl;
		if (x != null)
		{
		  for (int j = 0; j < DIM_SPACE; j++)
		  {
			x[i * sx + j] = _x[j];
		  }
		}
		if (dxdxi != null)
		{
		  for (int j = 0; j < DIM_SPACE; j++)
		  {
			for (int k = 0; k < 2; k++)
			{
			  dxdxi[i * sdxdxi + 2 * j + k] = _dxdxi(j,k);
			}
		  }
		}
	}
	if (ok)
	{
		return;
	}


// THESE LAST LINES ARE COPIED FROM CurvedElements::CalcSurfaceTransformation

	ArrayMem<Vec<DIM_SPACE>,100> coefs = new ArrayMem<Vec<DIM_SPACE>,100>(info.ndof);
	GetCoefficients(info, coefs);

	ArrayMem<T, 100> shapes_mem = new ArrayMem<T, 100>(info.ndof);
	TFlatVector<T> shapes = new TFlatVector<T>(info.ndof, shapes_mem[0]);

	ArrayMem<T, 100> dshapes_mem = new ArrayMem<T, 100>(info.ndof * 2);
	MatrixFixWidth<2,T> dshapes = new MatrixFixWidth<2,T>(info.ndof, shapes_mem[0]);



	if (x != null)
	{
	if (info.order == 1 && type == TRIG)
	{
		for (int j = 0; j < npts; j++)
		{
		Point<2,T> vxi = new Point<2,T>(xi[j * sxi], xi[j * sxi + 1]);

		Point<DIM_SPACE,T> val = new Point<DIM_SPACE,T>();
				for (int k = 0; k < DIM_SPACE; k++)
				{
				  val(k) = coefs[2](k) + (coefs[0](k) - coefs[2](k)) * vxi(0) + (coefs[1](k) - coefs[2](k)) * vxi(1);
				}
				/*
				(coefs[2]);
		val += (coefs[0]-coefs[2]) * vxi(0);
		val += (coefs[1]-coefs[2]) * vxi(1);
				*/
		for (int k = 0; k < DIM_SPACE; k++)
		{
		  x[j * sx + k] = val(k);
		}
		}
	}
	else
	{
	  for (int j = 0; j < npts; j++)
	  {
		  Point<2,T> vxi = new Point<2,T>(xi[j * sxi], xi[j * sxi + 1]);
		  CalcElementShapes(info, vxi, shapes);

		  Point<DIM_SPACE,T> val = T(0.0);
		  for (int i = 0; i < coefs.Size(); i++)
		  {
				for (int k = 0; k < DIM_SPACE; k++)
				{
				  val(k) += shapes(i) * coefs[i](k);
				}
		  }

		  for (int k = 0; k < DIM_SPACE; k++)
		  {
		x[j * sx + k] = val(k);
		  }
	  }
	}
	}

	if (dxdxi != null)
	{
	if (info.order == 1 && type == TRIG)
	{
		Point<2,T> xij = new Point<2,T>(xi[0], xi[1]);
		CalcElementDShapes(info, xij, dshapes);

		Mat<3,2,T> dxdxij = new Mat<3,2,T>();
		dxdxij = 0.0;
		for (int i = 0; i < coefs.Size(); i++)
		{
		  for (int j = 0; j < DIM_SPACE; j++)
		  {
		for (int k = 0; k < 2; k++)
		{
		  dxdxij(j,k) += dshapes(i,k) * coefs[i](j);
		}
		  }
		}


		for (int ip = 0; ip < npts; ip++)
		{
		  for (int j = 0; j < DIM_SPACE; j++)
		  {
		for (int k = 0; k < 2; k++)
		{
		  dxdxi[ip * sdxdxi + 2 * j + k] = dxdxij(j,k);
		}
		  }
		}
	}
	else
	{
		for (int j = 0; j < npts; j++)
		{
		Point<2,T> vxi = new Point<2,T>(xi[j * sxi], xi[j * sxi + 1]);
		CalcElementDShapes(info, vxi, dshapes);

		Mat<DIM_SPACE,2,T> ds = new Mat<DIM_SPACE,2,T>();
		ds = 0.0;
		for (int i = 0; i < coefs.Size(); i++)
		{
		  for (int j = 0; j < DIM_SPACE; j++)
		  {
			for (int k = 0; k < 2; k++)
			{
			  ds(j,k) += dshapes(i,k) * coefs[i](j);
			}
		  }
		}
		// (*dxdxi)[ip] = ds;

		for (int k = 0; k < 2 * DIM_SPACE; k++)
		{
		  dxdxi[j * sdxdxi + k] = ds(k);
		}
		}
	}
	}
  }

  public void CalcMultiPointElementTransformation(Array< Point < 3>> xi, ElementIndex elnr, Array< Point < 3>> x, Array< Mat < 3,3>> dxdxi)
  {
//C++ TO C# CONVERTER TODO TASK: C# does not have an equivalent to pointers to value types:
//ORIGINAL LINE: double * px = (x) ? &(*x)[0](0) : null;
	double px = (x) != null ?  x[0](0) : null;
//C++ TO C# CONVERTER TODO TASK: C# does not have an equivalent to pointers to value types:
//ORIGINAL LINE: double * pdxdxi = (dxdxi) ? &(*dxdxi)[0](0) : null;
	double pdxdxi = (dxdxi) != null ?  dxdxi[0](0) : null;

	CalcMultiPointElementTransformation(elnr, xi.Size(), xi[0](0), 3, px, 3, pdxdxi, 9);

	return;
#if OLD

	if (mesh.coarsemesh)
	{
	HPRefElement hpref_el = (*mesh.hpelements)[mesh[elnr].hp_elnr];

	// xi umrechnen
	double[] lami = new double[8];
	FlatVector vlami = new FlatVector(8, ref lami);


	ArrayMem<Point < 3>, 50> coarse_xi = new ArrayMem<Point < 3>, 50>(xi.Size());

	for (int pi = 0; pi < xi.Size(); pi++)
	{
		vlami = 0;
		mesh[elnr].GetShapeNew(xi[pi], vlami);

		Point < 3> cxi(0,0,0);
		for (int i = 0; i < hpref_el.np; i++)
		{
		  for (int j = 0; j < 3; j++)
		  {
		cxi(j) += hpref_el.param[i][j] * lami[i];
		  }
		}

		coarse_xi[pi] = cxi;
	}

	mesh.coarsemesh.GetCurvedElements().CalcMultiPointElementTransformation(coarse_xi, hpref_el.coarse_elnr, x, dxdxi);


	Mat < 3,3> trans, dxdxic;
	if (dxdxi != null)
	{
		MatrixFixWidth < 3> dlami(8);
		dlami = 0;

		for (int pi = 0; pi < xi.Size(); pi++)
		{
		mesh[elnr].GetDShapeNew(xi[pi], dlami);

		trans = 0;
		for (int k = 0; k < 3; k++)
		{
		  for (int l = 0; l < 3; l++)
		  {
			for (int i = 0; i < hpref_el.np; i++)
			{
			  trans(l,k) += hpref_el.param[i][l] * dlami(i, k);
			}
		  }
		}

		dxdxic = dxdxi[pi];
		dxdxi[pi] = dxdxic * trans;
		}
	}

	return;
	}








	Vector shapes = new Vector();
	MatrixFixWidth < 3> dshapes;


	Element el = mesh[elnr];
	ELEMENT_TYPE type = el.GetType();

	ElementInfo info = new ElementInfo();
	info.elnr = elnr;
	info.order = order;
	info.ndof = info.nv = MeshTopology.GetNPoints(type);
	if (info.order > 1)
	{
	MeshTopology top = mesh.GetTopology();

	info.nedges = top.GetElementEdges(elnr + 1, info.edgenrs, 0);
	for (int i = 0; i < info.nedges; i++)
	{
	  info.edgenrs[i]--;
	}

	info.nfaces = top.GetElementFaces(elnr + 1, info.facenrs, 0);
	for (int i = 0; i < info.nfaces; i++)
	{
	  info.facenrs[i]--;
	}

	for (int i = 0; i < info.nedges; i++)
	{
	  info.ndof += edgecoeffsindex[info.edgenrs[i] + 1] - edgecoeffsindex[info.edgenrs[i]];
	}
	for (int i = 0; i < info.nfaces; i++)
	{
	  info.ndof += facecoeffsindex[info.facenrs[i] + 1] - facecoeffsindex[info.facenrs[i]];
	}
	// info.ndof += facecoeffsindex[info.facenr+1] - facecoeffsindex[info.facenr];
	}

	Array<Vec < 3>> coefs = new Array<Vec < 3>>(info.ndof);
	GetCoefficients(info, coefs[0]);
	if (x != null)
	{
	for (int j = 0; j < xi.Size(); j++)
	{
		CalcElementShapes(info, xi[j], shapes);
		x[j] = 0;
		for (int i = 0; i < coefs.Size(); i++)
		{
		  x[j] += shapes(i) * coefs[i];
		}
	}
	}

	if (dxdxi != null)
	{
	if (info.order == 1 && type == TET)
	{
		if (xi.Size() > 0)
		{
		CalcElementDShapes(info, xi[0], dshapes);
		Mat < 3,3> ds;
		ds = 0;
		for (int i = 0; i < coefs.Size(); i++)
		{
		  for (int j = 0; j < 3; j++)
		  {
			for (int k = 0; k < 3; k++)
			{
			  ds(j,k) += dshapes(i,k) * coefs[i](j);
			}
		  }
		}

		for (int ip = 0; ip < xi.Size(); ip++)
		{
		  dxdxi[ip] = ds;
		}
		}
	}
	else
	{
	  for (int ip = 0; ip < xi.Size(); ip++)
	  {
		  CalcElementDShapes(info, xi[ip], dshapes);

		  Mat < 3,3> ds;
		  ds = 0;
		  for (int i = 0; i < coefs.Size(); i++)
		  {
		for (int j = 0; j < 3; j++)
		{
		  for (int k = 0; k < 3; k++)
		  {
			ds(j,k) += dshapes(i,k) * coefs[i](j);
		  }
		}
		  }
		  dxdxi[ip] = ds;
	  }
	}
	}
#endif
  }

//C++ TO C# CONVERTER TODO TASK: The original C++ template specifier was replaced with a C# generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <typename T>

  // extern int multipointtrafovar;
  public void CalcMultiPointElementTransformation<T>(ElementIndex elnr, int n, T[] xi, uint sxi, T[] x, uint sx, T[] dxdxi, uint sdxdxi)
  {
	// multipointtrafovar++;
	/*
	static int timer = NgProfiler::CreateTimer ("calcmultipointelementtrafo");
	static int timer1 = NgProfiler::CreateTimer ("calcmultipointelementtrafo 1");
	static int timer2 = NgProfiler::CreateTimer ("calcmultipointelementtrafo 2");
	static int timer3 = NgProfiler::CreateTimer ("calcmultipointelementtrafo 3");
	static int timer4 = NgProfiler::CreateTimer ("calcmultipointelementtrafo 4");
	static int timer5 = NgProfiler::CreateTimer ("calcmultipointelementtrafo 5");
	NgProfiler::RegionTimer reg(timer);
	*/
	// NgProfiler::StartTimer (timer);
	// NgProfiler::StartTimer (timer1);
	if (mesh.coarsemesh)
	{
	HPRefElement hpref_el = (*mesh.hpelements)[mesh[elnr].hp_elnr];

	// xi umrechnen
	T[] lami = Arrays.InitializeWithDefaultInstances<T>(8);
	TFlatVector<T> vlami = new TFlatVector<T>(8, lami[0]);


	ArrayMem<T, 100> coarse_xi = new ArrayMem<T, 100>((uint)(3 * n));

	for (int pi = 0; pi < n; pi++)
	{
		vlami = 0;
		Point<3,T> pxi = new Point<3,T>();
		for (int j = 0; j < 3; j++)
		{
		  pxi(j) = xi[pi * sxi + j];
		}

		mesh[elnr].GetShapeNew(pxi, vlami);

		Point<3,T> cxi = new Point<3,T>(0, 0, 0);
		for (int i = 0; i < hpref_el.np; i++)
		{
		  for (int j = 0; j < 3; j++)
		  {
		cxi(j) += hpref_el.param[i][j] * lami[i];
		  }
		}

		for (int j = 0; j < 3; j++)
		{
		  coarse_xi[3 * pi + j] = cxi(j);
		}
	}

	mesh.coarsemesh.GetCurvedElements().CalcMultiPointElementTransformation(hpref_el.coarse_elnr, n, coarse_xi[0], 3, x, sx, dxdxi, sdxdxi);

	Mat<3,3,T> trans = new Mat<3,3,T>();
	Mat<3,3,T> dxdxic = new Mat<3,3,T>();
	if (dxdxi != null)
	{
		MatrixFixWidth<3,T> dlami = new MatrixFixWidth<3,T>(8);
		dlami = T(0);

		for (int pi = 0; pi < n; pi++)
		{
		Point<3,T> pxi = new Point<3,T>();
		for (int j = 0; j < 3; j++)
		{
		  pxi(j) = xi[pi * sxi + j];
		}

		mesh[elnr].GetDShapeNew(pxi, dlami);

		trans = 0;
		for (int k = 0; k < 3; k++)
		{
		  for (int l = 0; l < 3; l++)
		  {
			for (int i = 0; i < hpref_el.np; i++)
			{
			  trans(l,k) += hpref_el.param[i][l] * dlami(i, k);
			}
		  }
		}

		Mat<3,3,T> mat_dxdxic = new Mat<3,3,T>();
		Mat<3,3,T> mat_dxdxi = new Mat<3,3,T>();
		for (int j = 0; j < 3; j++)
		{
		  for (int k = 0; k < 3; k++)
		  {
			mat_dxdxic(j,k) = dxdxi[pi * sdxdxi + 3 * j + k];
		  }
		}

		mat_dxdxi = mat_dxdxic * trans;

		for (int j = 0; j < 3; j++)
		{
		  for (int k = 0; k < 3; k++)
		  {
			dxdxi[pi * sdxdxi + 3 * j + k] = mat_dxdxi(j,k);
		  }
		}

		// dxdxic = (*dxdxi)[pi];
		// (*dxdxi)[pi] = dxdxic * trans;
		}
	}
	return;
	}

	// NgProfiler::StopTimer (timer1);
	// NgProfiler::StartTimer (timer2);


	Element el = mesh[elnr];
	ELEMENT_TYPE type = el.GetType();


	ElementInfo info = new ElementInfo();
	info.elnr = elnr;
	info.order = order;
	info.ndof = info.nv = MeshTopology.GetNPoints(type);
	if (info.order > 1)
	{
	MeshTopology top = mesh.GetTopology();

	info.nedges = top.GetElementEdges(elnr + 1, info.edgenrs, 0);
	for (int i = 0; i < info.nedges; i++)
	{
	  info.edgenrs[i]--;
	}

	info.nfaces = top.GetElementFaces(elnr + 1, info.facenrs, 0);
	for (int i = 0; i < info.nfaces; i++)
	{
	  info.facenrs[i]--;
	}

	for (int i = 0; i < info.nedges; i++)
	{
	  info.ndof += edgecoeffsindex[info.edgenrs[i] + 1] - edgecoeffsindex[info.edgenrs[i]];
	}
	for (int i = 0; i < info.nfaces; i++)
	{
	  info.ndof += facecoeffsindex[info.facenrs[i] + 1] - facecoeffsindex[info.facenrs[i]];
	}
	// info.ndof += facecoeffsindex[info.facenr+1] - facecoeffsindex[info.facenr];
	}

	// NgProfiler::StopTimer (timer2);
	// NgProfiler::StartTimer (timer3);


	bool ok = true;
	for (int i = 0; i < n; i++)
	{
		Point<3,T> _xi = new Point<3,T>(xi[i * sxi], xi[i * sxi + 1], xi[i * sxi + 2]);
		Point<3,T> _x = new Point<3,T>();
		Mat<3,3,T> _dxdxi = new Mat<3,3,T>();
		if (!EvaluateMapping(info, _xi, _x, _dxdxi))
		{
			  ok = false;
			  break;
		}
		// cout << "x = " << _x << ", dxdxi = " << _dxdxi << endl;
		if (x != null)
		{
		  for (int j = 0; j < 3; j++)
		  {
			x[i * sx + j] = _x[j];
		  }
		}
		if (dxdxi != null)
		{
		  for (int j = 0; j < 3; j++)
		  {
			for (int k = 0; k < 3; k++)
			{
			  dxdxi[i * sdxdxi + 3 * j + k] = _dxdxi(j,k);
			}
		  }
		}
	}
	if (ok)
	{
		return;
	}

	ArrayMem<Vec < 3>,100> coefs = new ArrayMem<Vec < 3>,100>(info.ndof);
	ArrayMem<T,500> shapes_mem = new ArrayMem<T,500>(info.ndof);

	TFlatVector<T> shapes = new TFlatVector<T>(info.ndof, shapes_mem[0]);

	ArrayMem<T,1500> dshapes_mem = new ArrayMem<T,1500>(3 * info.ndof);
	MatrixFixWidth<3,T> dshapes = new MatrixFixWidth<3,T>(info.ndof, dshapes_mem[0]);

	// NgProfiler::StopTimer (timer3);
	// NgProfiler::StartTimer (timer4);

	GetCoefficients(info, coefs[0]);
	if (x != null)
	{
	for (int j = 0; j < n; j++)
	{
		Point<3,T> xij = new Point<3,T>();
		Point<3,T> xj = new Point<3,T>();
		for (int k = 0; k < 3; k++)
		{
		  xij(k) = xi[j * sxi + k];
		}
		CalcElementShapes(info, xij, shapes);
		xj = T(0.0);
		for (int i = 0; i < coefs.Size(); i++)
		{
			  for (int k = 0; k < 3; k++)
			  {
				xj(k) += shapes(i) * coefs[i](k);
			  }
		}

			// cout << "old, xj = " << xj << endl;

		for (int k = 0; k < 3; k++)
		{
		  x[j * sx + k] = xj(k);
		}
	}
	}


	// NgProfiler::StopTimer (timer4);
	// NgProfiler::StartTimer (timer5);

	if (dxdxi != null)
	{
	if (info.order == 1 && type == TET)
	{
		if (n > 0)
		{

		Point<3,T> xij = new Point<3,T>();
		for (int k = 0; k < 3; k++)
		{
		  xij(k) = xi[k];
		}

		CalcElementDShapes(info, xij, dshapes);

		Mat<3,3,T> dxdxij = new Mat<3,3,T>();
		dxdxij = 0.0;
		for (int i = 0; i < coefs.Size(); i++)
		{
		  for (int j = 0; j < 3; j++)
		  {
			for (int k = 0; k < 3; k++)
			{
			  dxdxij(j,k) += dshapes(i,k) * coefs[i](j);
			}
		  }
		}


		for (int ip = 0; ip < n; ip++)
		{
		  for (int j = 0; j < 3; j++)
		  {
			for (int k = 0; k < 3; k++)
			{
			  dxdxi[ip * sdxdxi + 3 * j + k] = dxdxij(j,k);
			}
		  }
		}
		}
	}
	else
	{
		for (int ip = 0; ip < n; ip++)
		{
		Point<3,T> xij = new Point<3,T>();
		for (int k = 0; k < 3; k++)
		{
		  xij(k) = xi[ip * sxi + k];
		}

				CalcElementDShapes(info, xij, dshapes);


		Mat<3,3,T> dxdxij = new Mat<3,3,T>();
		dxdxij = 0.0;
		for (int i = 0; i < coefs.Size(); i++)
		{
		  for (int j = 0; j < 3; j++)
		  {
			for (int k = 0; k < 3; k++)
			{
			  dxdxij(j,k) += dshapes(i,k) * coefs[i](j);
			}
		  }
		}

				// cout << "old, jac = " << dxdxij << endl;

		for (int j = 0; j < 3; j++)
		{
		  for (int k = 0; k < 3; k++)
		  {
			dxdxi[ip * sdxdxi + 3 * j + k] = dxdxij(j,k);
		  }
		}

				/*
				T dxdxi00 = T(0.0);
				T dxdxi01 = T(0.0);
				T dxdxi02 = T(0.0);
				T dxdxi10 = T(0.0);
				T dxdxi11 = T(0.0);
				T dxdxi12 = T(0.0);
				T dxdxi20 = T(0.0);
				T dxdxi21 = T(0.0);
				T dxdxi22 = T(0.0);
				
		for (int i = 0; i < coefs.Size(); i++)
				  {
				    T ds0 = dshapes(i,0);
				    T ds1 = dshapes(i,1);
				    T ds2 = dshapes(i,2);
				    T cf0 = coefs[i](0);
				    T cf1 = coefs[i](1);
				    T cf2 = coefs[i](2);

				    dxdxi00 += ds0*cf0;
				    dxdxi01 += ds1*cf0;
				    dxdxi02 += ds2*cf0;
				    dxdxi10 += ds0*cf1;
				    dxdxi11 += ds1*cf1;
				    dxdxi12 += ds2*cf1;
				    dxdxi20 += ds0*cf2;
				    dxdxi21 += ds1*cf2;
				    dxdxi22 += ds2*cf2;
				  }

				dxdxi[ip*sdxdxi+3*0+0] = dxdxi00;
				dxdxi[ip*sdxdxi+3*0+1] = dxdxi01;
				dxdxi[ip*sdxdxi+3*0+2] = dxdxi02;

				dxdxi[ip*sdxdxi+3*1+0] = dxdxi10;
				dxdxi[ip*sdxdxi+3*1+1] = dxdxi11;
				dxdxi[ip*sdxdxi+3*1+2] = dxdxi12;
				
				dxdxi[ip*sdxdxi+3*2+0] = dxdxi20;
				dxdxi[ip*sdxdxi+3*2+1] = dxdxi21;
				dxdxi[ip*sdxdxi+3*2+2] = dxdxi22;
				*/
		}
	}
	}
	// NgProfiler::StopTimer (timer5);
	// NgProfiler::StopTimer (timer);
  }





//C++ TO C# CONVERTER TODO TASK: The original C++ template specifier was replaced with a C# generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <typename T>
  private void CalcSegmentTransformation<T>(T xi, SegmentIndex elnr, Point<3,T> x = null, Vec<3,T> dxdxi = null)
  {
	  CalcSegmentTransformation(xi, elnr, null, null, null);
  }
  private void CalcSegmentTransformation<T>(T xi, SegmentIndex elnr, Point<3,T> x = null)
  {
	  CalcSegmentTransformation(xi, elnr, null, null, null);
  }
  private void CalcSegmentTransformation<T>(T xi, SegmentIndex elnr)
  {
	  CalcSegmentTransformation(xi, elnr, null, null, null);
  }
  private void CalcSegmentTransformation<T>(T xi, SegmentIndex elnr, Point<3,T> x = null, Vec<3,T> dxdxi = null, ref bool curved)
  {
	if (mesh.coarsemesh)
	{
	HPRefElement hpref_el = (*mesh.hpelements)[mesh[elnr].hp_elnr];

	// xi umrechnen
	T[] lami = {xi, 1 - xi};
	T[] dlami = {1, -1};

	T coarse_xi = 0;
	T trans = 0;
	for (int i = 0; i < 2; i++)
	{
		coarse_xi += hpref_el.param[i][0] * lami[i];
		trans += hpref_el.param[i][0] * dlami[i];
	}

	mesh.coarsemesh.GetCurvedElements().CalcSegmentTransformation(coarse_xi, hpref_el.coarse_elnr, x, dxdxi, curved);
	if (dxdxi != null)
	{
		dxdxi *= trans;
	}

	return;
	}



	// TVector<T> shapes, dshapes;
	//     Array<Vec<3> > coefs;

	SegmentInfo info = new SegmentInfo();
	info.elnr = elnr;
	info.order = order;
	info.ndof = info.nv = 2;

	if (info.order > 1)
	{
	MeshTopology top = mesh.GetTopology();
	info.edgenr = top.GetSegmentEdge(elnr + 1) - 1;
	info.ndof += edgeorder[info.edgenr] - 1;
	}

	ArrayMem<Vec < 3>,100> coefs = new ArrayMem<Vec < 3>,100>(info.ndof);
	ArrayMem<T, 100> shapes_mem = new ArrayMem<T, 100>(info.ndof);
	TFlatVector<T> shapes = new TFlatVector<T>(info.ndof, shapes_mem[0]);
	ArrayMem<T, 200> dshapes_mem = new ArrayMem<T, 200>(info.ndof);
	TFlatVector<T> dshapes = new TFlatVector<T>(info.ndof, dshapes_mem[0]);


	CalcElementShapes(info, xi, shapes);
	GetCoefficients(info, coefs);

	x = null;
	for (int i = 0; i < shapes.Size(); i++)
	{
	  // *x += shapes(i) * coefs[i];
	  for (int j = 0; j < 3; j++)
	  {
		x(j) += shapes(i) * coefs[i](j);
	  }
	}


	if (dxdxi != null)
	{
	CalcElementDShapes(info, xi, dshapes);

	dxdxi = null;
	for (int i = 0; i < shapes.Size(); i++)
	{
	  for (int j = 0; j < 3; j++)
	  {
		dxdxi(j) += dshapes(i) * coefs[i](j);
	  }
	}
	}

	if (curved)
	{
	  curved = (info.order > 1);
	}

	// cout << "Segment, |x| = " << Abs2(Vec<3> (*x) ) << endl;
  }

  private void CalcSurfaceTransformation(Point < 2> xi, SurfaceElementIndex elnr, Point < 3> * x = null, Mat<3, 2> dxdxi = null)
  {
	  CalcSurfaceTransformation(xi, elnr, null, null, null);
  }
  private void CalcSurfaceTransformation(Point < 2> xi, SurfaceElementIndex elnr, Point < 3> * x = null)
  {
	  CalcSurfaceTransformation(xi, elnr, null, null, null);
  }
  private void CalcSurfaceTransformation(Point < 2> xi, SurfaceElementIndex elnr)
  {
	  CalcSurfaceTransformation(xi, elnr, null, null, null);
  }
//C++ TO C# CONVERTER TODO TASK: Pointer arithmetic is detected on the parameter 'x', so pointers on this parameter are left unchanged:
  private void CalcSurfaceTransformation(Point < 2> xi, SurfaceElementIndex elnr, Point < 3> * x = null, Mat<3, 2> dxdxi = null, ref bool curved)
  {
	if (mesh.coarsemesh)
	{
	HPRefElement hpref_el = (*mesh.hpelements)[mesh[elnr].hp_elnr];

	// xi umrechnen
	double[] lami = new double[4];
	FlatVector vlami = new FlatVector(4, ref lami);
	vlami = 0;
	mesh[elnr].GetShapeNew(xi, vlami);

	Mat < 2,2> trans;
	Mat < 3,2> dxdxic;
	if (dxdxi != null)
	{
		MatrixFixWidth < 2> dlami(4);
		dlami = 0;
		mesh[elnr].GetDShapeNew(xi, dlami);

		trans = 0;
		for (int k = 0; k < 2; k++)
		{
		  for (int l = 0; l < 2; l++)
		  {
		for (int i = 0; i < hpref_el.np; i++)
		{
		  trans(l,k) += hpref_el.param[i][l] * dlami(i, k);
		}
		  }
		}
	}

	Point < 2> coarse_xi(0,0);
	for (int i = 0; i < hpref_el.np; i++)
	{
	  for (int j = 0; j < 2; j++)
	  {
		coarse_xi(j) += hpref_el.param[i][j] * lami[i];
	  }
	}

	mesh.coarsemesh.GetCurvedElements().CalcSurfaceTransformation(coarse_xi, hpref_el.coarse_elnr, x, dxdxic, curved);

	if (dxdxi != null)
	{
	  dxdxi = dxdxic * trans;
	}

	return;
	}




	Element2d el = mesh[elnr];
	ELEMENT_TYPE type = el.GetType();

	SurfaceElementInfo info = new SurfaceElementInfo();
	info.elnr = elnr;
	info.order = order;

	switch (type)
	{
	  case TRIG :
		  info.nv = 3;
		  break;
	  case QUAD :
		  info.nv = 4;
		  break;
	  case TRIG6:
		  info.nv = 6;
		  break;
	  case QUAD8 :
		  info.nv = 8;
		  break;
	  default:
	cerr << "undef element in CalcSurfaceTrafo" << "\n";
	break;
	}
	info.ndof = info.nv;

	if (info.order > 1)
	{
	MeshTopology top = mesh.GetTopology();

	top.GetSurfaceElementEdges(elnr + 1, info.edgenrs);
	for (int i = 0; i < info.edgenrs.Size(); i++)
	{
	  info.edgenrs[i]--;
	}
	info.facenr = top.GetSurfaceElementFace(elnr + 1) - 1;


	bool firsttry = true;
	bool problem = false;

	while (firsttry || problem)
	{
		problem = false;

		for (int i = 0; !problem && i < info.edgenrs.Size(); i++)
		{
		if (info.edgenrs[i] + 1 >= edgecoeffsindex.Size())
		{
		  problem = true;
		}
		else
		{
		  info.ndof += edgecoeffsindex[info.edgenrs[i] + 1] - edgecoeffsindex[info.edgenrs[i]];
		}
		}
		if (info.facenr + 1 >= facecoeffsindex.Size())
		{
		  problem = true;
		}
		else
		{
		  info.ndof += facecoeffsindex[info.facenr + 1] - facecoeffsindex[info.facenr];
		}

		if (problem && !firsttry)
		{
		  throw new Exception("something wrong with curved elements");
		}

		if (problem)
		{
		  BuildCurvedElements(null,order,rational);
		}

		firsttry = false;
	}
	}


	Point < 2> _xi(xi);
	Point < 3> _x;
	Mat < 3,2> _dxdxi;
	if (EvaluateMapping(info, _xi, _x, _dxdxi))
	{
		if (x != null)
		{
			*x = _x;
		}
		if (dxdxi != null)
		{
			dxdxi = _dxdxi;
		}
		return;
	}


	ArrayMem<Vec < 3>,100> coefs = new ArrayMem<Vec < 3>,100>(info.ndof);
	ArrayMem<double, 100> shapes_mem = new ArrayMem<double, 100>(info.ndof);
	TFlatVector<double> shapes = new TFlatVector<double>(info.ndof, shapes_mem[0]);
	ArrayMem<double, 200> dshapes_mem = new ArrayMem<double, 200>(2 * info.ndof);
	MatrixFixWidth < 2> dshapes(info.ndof, dshapes_mem[0]);


	CalcElementShapes(info, xi, shapes);
	GetCoefficients(info, coefs);

	*x = 0;
	for (int i = 0; i < coefs.Size(); i++)
	{
	  *x += shapes(i) * coefs[i];
	}

	if (dxdxi != null)
	{
	CalcElementDShapes(info, xi, dshapes);

	dxdxi = null;
	for (int i = 0; i < coefs.Size(); i++)
	{
	  for (int j = 0; j < 3; j++)
	  {
		for (int k = 0; k < 2; k++)
		{
		  dxdxi(j,k) += dshapes(i,k) * coefs[i](j);
		}
	  }
	}
	}

	if (curved)
	{
	  curved = (info.ndof > info.nv);
	}
  }

//C++ TO C# CONVERTER TODO TASK: Pointer arithmetic is detected on the parameter 'x', so pointers on this parameter are left unchanged:
  private void CalcElementTransformation(Point < 3> xi, ElementIndex elnr, Point < 3> * x = null, Mat<3, 3> dxdxi = null, object buffer = null, bool valid = false)
  {
	if (mesh.coarsemesh)
	{
	HPRefElement hpref_el = (*mesh.hpelements)[mesh[elnr].hp_elnr];

	// xi umrechnen
	double[] lami = new double[8];
	FlatVector vlami = new FlatVector(8, ref lami);
	vlami = 0;
	mesh[elnr].GetShapeNew<double> (xi, vlami);

	Mat < 3,3> trans, dxdxic;
	if (dxdxi != null)
	{
		MatrixFixWidth < 3> dlami(8);
		dlami = 0;
		mesh[elnr].GetDShapeNew(xi, dlami);

		trans = 0;
		for (int k = 0; k < 3; k++)
		{
		  for (int l = 0; l < 3; l++)
		  {
		for (int i = 0; i < hpref_el.np; i++)
		{
		  trans(l,k) += hpref_el.param[i][l] * dlami(i, k);
		}
		  }
		}
	}

	Point < 3> coarse_xi(0,0,0);
	for (int i = 0; i < hpref_el.np; i++)
	{
	  for (int j = 0; j < 3; j++)
	  {
		coarse_xi(j) += hpref_el.param[i][j] * lami[i];
	  }
	}

	mesh.coarsemesh.GetCurvedElements().CalcElementTransformation(coarse_xi, hpref_el.coarse_elnr, x, dxdxic);

	if (dxdxi != null)
	{
	  dxdxi = dxdxic * trans;
	}

	return;
	}


	Element el = mesh[elnr];
	ELEMENT_TYPE type = el.GetType();

	ElementInfo hinfo = new ElementInfo();
	ElementInfo info = (buffer) ? *(ElementInfo) buffer : hinfo;


	if (!valid)
	{
	info.elnr = elnr;
	info.order = order;
	info.ndof = info.nv = MeshTopology.GetNPoints(type);
	if (info.order > 1)
	{
		MeshTopology top = mesh.GetTopology();

		info.nedges = top.GetElementEdges(elnr + 1, info.edgenrs, 0);
		for (int i = 0; i < info.nedges; i++)
		{
		  info.edgenrs[i]--;
		}

		info.nfaces = top.GetElementFaces(elnr + 1, info.facenrs, 0);
		for (int i = 0; i < info.nfaces; i++)
		{
		  info.facenrs[i]--;
		}

		for (int i = 0; i < info.nedges; i++)
		{
		  info.ndof += edgecoeffsindex[info.edgenrs[i] + 1] - edgecoeffsindex[info.edgenrs[i]];
		}
		for (int i = 0; i < info.nfaces; i++)
		{
		  info.ndof += facecoeffsindex[info.facenrs[i] + 1] - facecoeffsindex[info.facenrs[i]];
		}
	}
	}

	ArrayMem<double,100> mem = new ArrayMem<double,100>(info.ndof);
	TFlatVector<double> shapes = new TFlatVector<double>(info.ndof, mem[0]);
	ArrayMem<double,100> dshapes_mem = new ArrayMem<double,100>(info.ndof * 3);
	MatrixFixWidth < 3> dshapes(info.ndof, dshapes_mem[0]);

	CalcElementShapes(info, xi, shapes);

	Vec < 3> * coefs = (info.ndof <= 10) ? info.hcoefs[0] : new Vec < 3> [info.ndof];

	if (info.ndof > 10 || !valid)
	{
	  GetCoefficients(info, coefs);
	}

	if (x != null)
	{
	*x = 0;
	for (int i = 0; i < shapes.Size(); i++)
	{
	  *x += shapes(i) * coefs[i];
	}
	}

	if (dxdxi != null)
	{
	if (valid && info.order == 1 && info.nv == 4) // a linear tet
	{
		dxdxi = info.hdxdxi;
	}
	else
	{
		CalcElementDShapes(info, xi, dshapes);

		dxdxi = null;
		for (int i = 0; i < shapes.Size(); i++)
		{
		  for (int j = 0; j < 3; j++)
		  {
		for (int k = 0; k < 3; k++)
		{
		  dxdxi(j,k) += dshapes(i,k) * coefs[i](j);
		}
		  }
		}

		info.hdxdxi = dxdxi;
	}
	}

	// *testout << "curved_elements, dshapes = " << endl << dshapes << endl;

	//    if (curved) *curved = (info.ndof > info.nv);

	if (info.ndof > 10)
	{
		coefs = null;
	}
  }






  private class SegmentInfo
  {
	public SegmentIndex elnr = new SegmentIndex();
	public int order;
	public int nv;
	public int ndof;
	public int edgenr;
  }

//C++ TO C# CONVERTER TODO TASK: The original C++ template specifier was replaced with a C# generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <typename T>
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void CalcElementShapes(SegmentInfo & info, T xi, TFlatVector<T> shapes) const
  private void CalcElementShapes<T>(SegmentInfo info, T xi, TFlatVector<T> shapes)
  {
	/*
	if (rational && info.order == 2)
	  {
	shapes.SetSize(3);
	double w = edgeweight[info.edgenr];
	shapes(0) = xi*xi;
	shapes(1) = (1-xi)*(1-xi);
	shapes(2) = 2*w*xi*(1-xi);
	shapes *= 1.0 / (1 + (w-1) *2*xi*(1-xi));
	return;
	  }
	*/

	// shapes.SetSize(info.ndof);
	shapes(0) = xi;
	shapes(1) = 1 - xi;

	if (info.order >= 2)
	{
	if (mesh[info.elnr][0] > mesh[info.elnr][1])
	{
	  xi = 1 - xi;
	}
	netgen.GlobalMembers.CalcEdgeShape(edgeorder[info.edgenr], 2 * xi - 1, shapes(2));
	}
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetCoefficients(SegmentInfo & info, Array<Vec<3>> & coefs) const
  private void GetCoefficients(SegmentInfo info, Array<Vec < 3>> coefs)
  {
	Segment el = mesh[info.elnr];

	coefs.SetSize(info.ndof);

	coefs[0] = Vec < 3> (mesh[el[0]]);
	coefs[1] = Vec < 3> (mesh[el[1]]);

	if (info.order >= 2)
	{
	int first = edgecoeffsindex[info.edgenr];
	int next = edgecoeffsindex[info.edgenr + 1];
	for (int i = 0; i < next - first; i++)
	{
	  coefs[i + 2] = edgecoeffs[first + i];
	}
	}
  }

//C++ TO C# CONVERTER TODO TASK: The original C++ template specifier was replaced with a C# generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <typename T>
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void CalcElementDShapes(SegmentInfo & info, T xi, TFlatVector<T> dshapes) const
  private void CalcElementDShapes<T>(SegmentInfo info, T xi, TFlatVector<T> dshapes)
  {
	/*
	if (rational && info.order == 2)
	  {
	dshapes.SetSize(3);
	double wi = edgeweight[info.edgenr];
	double shapes[3];
	shapes[0] = xi*xi;
	shapes[1] = (1-xi)*(1-xi);
	shapes[2] = 2*wi*xi*(1-xi);
	double w = 1 + (wi-1) *2*xi*(1-xi);
	double dw = (wi-1) * (2 - 4*xi);
	
	dshapes(0) = 2*xi;
	dshapes(1) = 2*(xi-1);
	dshapes(2) = 2*wi*(1-2*xi);

	for (int j = 0;j < 3; j++)
	  dshapes(j) = dshapes(j) / w - shapes[j] * dw / (w*w);
	return;
	  }
	*/



	// dshapes.SetSize(info.ndof);
	dshapes = 0;
	dshapes(0) = 1;
	dshapes(1) = -1;

	// int order = edgeorder[info.edgenr];

	if (info.order >= 2)
	{
	T fac = 2;
	if (mesh[info.elnr][0] > mesh[info.elnr][1])
	{
		xi = 1 - xi;
		fac *= -1;
	}
	netgen.GlobalMembers.CalcEdgeDx(edgeorder[info.edgenr], 2 * xi - 1, dshapes(2));
	for (int i = 2; i < dshapes.Size(); i++)
	{
	  dshapes(i) *= fac;
	}
	}

	// ??? not implemented ????
  }


  private class ElementInfo
  {
	public ElementIndex elnr = new ElementIndex();
	public int order;
	public int nv;
	public int ndof;
	public int nedges;
	public int nfaces;
	public int[] edgenrs = new int[12];
	public int[] facenrs = new int[6];
	public Mat < 3> hdxdxi;
	public Vec[] < 3> hcoefs = Arrays.InitializeWithDefaultInstances<Vec>(10); // enough for second order tets
  }

//C++ TO C# CONVERTER TODO TASK: The original C++ template specifier was replaced with a C# generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <typename T>
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void CalcElementShapes(ElementInfo & info, Point<3,T> xi, TFlatVector<T> shapes) const
  private void CalcElementShapes<T>(ElementInfo info, Point<3,T> xi, TFlatVector<T> shapes)
  {
	Element el = mesh[info.elnr];

	if (rational && info.order >= 2)
	{
	// shapes.SetSize(10);
	T w = 1;
	T[] lami = {xi(0), xi(1), xi(2), 1 - xi(0) - xi(1) - xi(2)};
	for (int j = 0; j < 4; j++)
	{
	  shapes(j) = lami[j] * lami[j];
	}

	ELEMENT_EDGE edges = MeshTopology.GetEdges1(TET);
	for (int j = 0; j < 6; j++)
	{
		double wi = edgeweight[info.edgenrs[j]];
		shapes(j + 4) = 2 * wi * lami[edges[j][0] - 1] * lami[edges[j][1] - 1];
		w += (wi - 1) * 2 * lami[edges[j][0] - 1] * lami[edges[j][1] - 1];
	}

	shapes *= 1.0 / w;
	return;
	}

	// shapes.SetSize(info.ndof);

	switch (el.GetType())
	{
	  case TET:
	  {
	  shapes(0) = xi(0);
	  shapes(1) = xi(1);
	  shapes(2) = xi(2);
	  shapes(3) = 1 - xi(0) - xi(1) - xi(2);

	  if (info.order == 1)
	  {
		  return;
	  }

	  int ii = 4;
	  ELEMENT_EDGE edges = MeshTopology.GetEdges1(TET);
	  for (int i = 0; i < 6; i++)
	  {
		  int eorder = edgeorder[info.edgenrs[i]];
		  if (eorder >= 2)
		  {
		  int vi1 = edges[i][0] - 1;
		  int vi2 = edges[i][1] - 1;
		  if (el[vi1] > el[vi2])
		  {
			  swap(vi1, vi2);
		  }

		  netgen.GlobalMembers.CalcScaledEdgeShape(eorder, shapes(vi1) - shapes(vi2), shapes(vi1) + shapes(vi2), shapes(ii));
		  ii += eorder - 1;
		  }
	  }
	  ELEMENT_FACE faces = MeshTopology.GetFaces1(TET);
	  for (int i = 0; i < 4; i++)
	  {
		  int forder = faceorder[info.facenrs[i]];
		  if (forder >= 3)
		  {
		  int[] fnums = {faces[i][0] - 1, faces[i][1] - 1, faces[i][2] - 1};
		  if (el[fnums[0]] > el[fnums[1]])
		  {
			  swap(fnums[0], fnums[1]);
		  }
		  if (el[fnums[1]] > el[fnums[2]])
		  {
			  swap(fnums[1], fnums[2]);
		  }
		  if (el[fnums[0]] > el[fnums[1]])
		  {
			  swap(fnums[0], fnums[1]);
		  }

		  netgen.GlobalMembers.CalcScaledTrigShape(forder, shapes(fnums[1]) - shapes(fnums[0]), shapes(fnums[2]), shapes(fnums[0]) + shapes(fnums[1]) + shapes(fnums[2]), shapes(ii));
		  ii += (forder - 1) * (forder - 2) / 2;
		  }
	  }

	  break;
	  }

	  case TET10:
	  {
	  T x = xi(0);
	  T y = xi(1);
	  T z = xi(2);
	  T lam4 = 1 - x - y - z;
	  /*
	    shapes(0) = xi(0);
	    shapes(1) = xi(1);
	    shapes(2) = xi(2);
	    shapes(3) = 1-xi(0)-xi(1)-xi(2);
	  */

	  shapes(0) = 2 * x * x - x;
	  shapes(1) = 2 * y * y - y;
	  shapes(2) = 2 * z * z - z;
	  shapes(3) = 2 * lam4 * lam4 - lam4;

	  shapes(4) = 4 * x * y;
	  shapes(5) = 4 * x * z;
	  shapes(6) = 4 * x * lam4;
	  shapes(7) = 4 * y * z;
	  shapes(8) = 4 * y * lam4;
	  shapes(9) = 4 * z * lam4;

	  break;
	  }

	  case PRISM:
	  {
	  T[] lami = {xi(0), xi(1), 1 - xi(0) - xi(1), xi(0), xi(1), 1 - xi(0) - xi(1)};
	  T[] lamiz = {1 - xi(2), 1 - xi(2), 1 - xi(2), xi(2), xi(2), xi(2)};
	  for (int i = 0; i < 6; i++)
	  {
		shapes(i) = lami[i] * lamiz[i];
	  }
	  for (int i = 6; i < info.ndof; i++)
	  {
		shapes(i) = 0;
	  }

	  if (info.order == 1)
	  {
		  return;
	  }


	  int ii = 6;
	  ELEMENT_EDGE edges = MeshTopology.GetEdges1(PRISM);
	  for (int i = 0; i < 6; i++) // horizontal edges
	  {
		  int eorder = edgeorder[info.edgenrs[i]];
		  if (eorder >= 2)
		  {
		  int vi1 = edges[i][0] - 1;
		  int vi2 = edges[i][1] - 1;
		  if (el[vi1] > el[vi2])
		  {
			  swap(vi1, vi2);
		  }

		  netgen.GlobalMembers.CalcScaledEdgeShape(eorder, lami[vi1] - lami[vi2], lami[vi1] + lami[vi2], shapes(ii));
		  T facz = (i < 3) ? (1 - xi(2)) : xi(2);
		  for (int j = 0; j < eorder - 1; j++)
		  {
			shapes(ii + j) *= facz;
		  }

		  ii += eorder - 1;
		  }
	  }

	  for (int i = 6; i < 9; i++) // vertical edges
	  {
		  int eorder = edgeorder[info.edgenrs[i]];
		  if (eorder >= 2)
		  {
		  int vi1 = edges[i][0] - 1;
		  int vi2 = edges[i][1] - 1;
		  if (el[vi1] > el[vi2])
		  {
			  swap(vi1, vi2);
		  }

		  T bubz = lamiz[vi1] * lamiz[vi2];
		  T polyz = lamiz[vi1] - lamiz[vi2];
		  T bubxy = lami[vi1];

		  for (int j = 0; j < eorder - 1; j++)
		  {
			  shapes(ii + j) = bubxy * bubz;
			  bubz *= polyz;
		  }
		  ii += eorder - 1;
		  }
	  }

	  // FACE SHAPES
	  ELEMENT_FACE faces = MeshTopology.GetFaces1(PRISM);
	  for (int i = 0; i < 2; i++)
	  {
		  int forder = faceorder[info.facenrs[i]];
		  if (forder < 3)
		  {
			  continue;
		  }
		  int[] fav = {faces[i][0] - 1, faces[i][1] - 1, faces[i][2] - 1};
		  if (el[fav[0]] > el[fav[1]])
		  {
			  swap(fav[0],fav[1]);
		  }
		  if (el[fav[1]] > el[fav[2]])
		  {
			  swap(fav[1],fav[2]);
		  }
		  if (el[fav[0]] > el[fav[1]])
		  {
			  swap(fav[0],fav[1]);
		  }

		  netgen.GlobalMembers.CalcTrigShape(forder, lami[fav[2]] - lami[fav[1]], lami[fav[0]], shapes(ii));

		  int ndf = (forder + 1) * (forder + 2) / 2 - 3 - 3 * (forder - 1);
		  for (int j = 0; j < ndf; j++)
		  {
		shapes(ii + j) *= lamiz[fav[1]];
		  }
		  ii += ndf;
	  }
	  break;
	  }

	  case PRISM15:
	  {
	  shapes = 0.0;
	  T x = xi(0);
	  T y = xi(1);
	  T z = xi(2);
		  T lam = 1 - x - y;
		  T lamz = 1 - z;
		  shapes[0] = (2 * x * x - x) * (2 * lamz * lamz - lamz);
		  shapes[1] = (2 * y * y - y) * (2 * lamz * lamz - lamz);
		  shapes[2] = (2 * lam * lam - lam) * (2 * lamz * lamz - lamz);
		  shapes[3] = (2 * x * x - x) * (2 * z * z - z);
		  shapes[4] = (2 * y * y - y) * (2 * z * z - z);
		  shapes[5] = (2 * lam * lam - lam) * (2 * z * z - z);
		  shapes[6] = 4 * x * y * (2 * lamz * lamz - lamz);
		  shapes[7] = 4 * x * lam * (2 * lamz * lamz - lamz);
		  shapes[8] = 4 * y * lam * (2 * lamz * lamz - lamz);
		  shapes[9] = x * 4 * z * (1 - z);
		  shapes[10] = y * 4 * z * (1 - z);
		  shapes[11] = lam * 4 * z * (1 - z);
		  shapes[12] = 4 * x * y * (2 * z * z - z);
		  shapes[13] = 4 * x * lam * (2 * z * z - z);
		  shapes[14] = 4 * y * lam * (2 * z * z - z);
		  break;
	  }

	  case PYRAMID:
	  {
	  shapes = 0.0;
	  T x = xi(0);
	  T y = xi(1);
	  T z = xi(2);

	  // if (z == 1.) z = 1-1e-10;
		  z *= (1 - 1e-12);
	  shapes[0] = (1 - z - x) * (1 - z - y) / (1 - z);
	  shapes[1] = x * (1 - z - y) / (1 - z);
	  shapes[2] = x * y / (1 - z);
	  shapes[3] = (1 - z - x) * y / (1 - z);
	  shapes[4] = z;

	  if (info.order == 1)
	  {
		  return;
	  }

		  T[] sigma = {sigma[0] = ((1 - z - x) + (1 - z - y)), sigma[1] = (x + (1 - z - y)), sigma[2] = (x + y), sigma[3] = ((1 - z - x) + y)};

	  int ii = 5;
	  ELEMENT_EDGE edges = MeshTopology.GetEdges1(PYRAMID);
	  for (int i = 0; i < 4; i++) // horizontal edges
	  {
		  int eorder = edgeorder[info.edgenrs[i]];
		  if (eorder >= 2)
		  {
		  int vi1 = (edges[i][0] - 1);
		  int vi2 = (edges[i][1] - 1);
		  if (el[vi1] > el[vi2])
		  {
			  swap(vi1, vi2);
		  }

				  netgen.GlobalMembers.CalcScaledEdgeShape(eorder, sigma[vi1] - sigma[vi2], 1 - z, shapes(ii));
		  T fac = (shapes[vi1] + shapes[vi2]) / (1 - z);
		  for (int j = 0; j < eorder - 1; j++)
		  {
			shapes(ii + j) *= fac;
		  }

		  ii += eorder - 1;
		  }
	  }



	  break;
	  }

	  case PYRAMID13:
	  {
	  shapes = 0.0;
	  T x = xi(0);
	  T y = xi(1);
	  T z = xi(2);
		  z *= 1 - 1e-12;
		  shapes[0] = (-z + z * (2 * x + z - 1) * (2 * y + z - 1) / (-z + 1) + (-2 * x - z + 2) * (-2 * y - z + 2)) * (-0.5 * x - 0.5 * y - 0.5 * z + 0.25);
		  shapes[1] = (0.5 * x - 0.5 * y - 0.25) * (-z - z * (2 * x + z - 1) * (2 * y + z - 1) / (-z + 1) + (2 * x + z) * (-2 * y - z + 2));
		  shapes[2] = (-z + z * (2 * x + z - 1) * (2 * y + z - 1) / (-z + 1) + (2 * x + z) * (2 * y + z)) * (0.5 * x + 0.5 * y + 0.5 * z - 0.75);
		  shapes[3] = (-0.5 * x + 0.5 * y - 0.25) * (-z - z * (2 * x + z - 1) * (2 * y + z - 1) / (-z + 1) + (2 * y + z) * (-2 * x - z + 2));
		  shapes[4] = z * (2 * z - 1);
		  shapes[5] = 2 * x * (-2 * x - 2 * z + 2) * (-2 * y - 2 * z + 2) / (-2 * z + 2);
		  shapes[6] = 4 * x * y * (-2 * x - 2 * z + 2) / (-2 * z + 2);
		  shapes[7] = 2 * y * (-2 * x - 2 * z + 2) * (-2 * y - 2 * z + 2) / (-2 * z + 2);
		  shapes[8] = 4 * x * y * (-2 * y - 2 * z + 2) / (-2 * z + 2);
		  shapes[9] = z * (-2 * x - 2 * z + 2) * (-2 * y - 2 * z + 2) / (-z + 1);
		  shapes[10] = 2 * x * z * (-2 * y - 2 * z + 2) / (-z + 1);
		  shapes[11] = 4 * x * y * z / (-z + 1);
		  shapes[12] = 2 * y * z * (-2 * x - 2 * z + 2) / (-z + 1);
		  break;
	  }

	  case HEX:
	  {
	  shapes = 0.0;
	  T x = xi(0);
	  T y = xi(1);
	  T z = xi(2);

	  shapes[0] = (1 - x) * (1 - y) * (1 - z);
	  shapes[1] = x * (1 - y) * (1 - z);
	  shapes[2] = x * y * (1 - z);
	  shapes[3] = (1 - x) * y * (1 - z);
	  shapes[4] = (1 - x) * (1 - y) * (z);
	  shapes[5] = x * (1 - y) * (z);
	  shapes[6] = x * y * (z);
	  shapes[7] = (1 - x) * y * (z);

	  if (info.order == 1)
	  {
		  return;
	  }

	  T[] mu = {(1 - x) + (1 - y) + (1 - z), x + (1 - y) + (1 - z), x + y + (1 - z), (1 - x) + y + (1 - z), (1 - x) + (1 - y) + (z), x + (1 - y) + (z), x + y + (z), (1 - x) + y + (z)};

	  int ii = 8;
	  ELEMENT_EDGE edges = MeshTopology.GetEdges1(HEX);

	  for (int i = 0; i < 8; i++)
	  {
		  int eorder = edgeorder[info.edgenrs[i]];
		  if (eorder >= 2)
		  {
		  int vi1 = edges[i][0] - 1;
		  int vi2 = edges[i][1] - 1;
		  if (el[vi1] > el[vi2])
		  {
			  swap(vi1, vi2);
		  }

		  netgen.GlobalMembers.CalcEdgeShape(eorder, mu[vi1] - mu[vi2], shapes(ii));
		  T lame = shapes(vi1) + shapes(vi2);
		  for (int j = 0; j < order - 1; j++)
		  {
			shapes(ii + j) *= lame;
		  }
		  ii += eorder - 1;
		  }
	  }


	  break;
	  }

	  case HEX20:
	  {
	  shapes = 0.0;
	  T x = xi(0);
	  T y = xi(1);
	  T z = xi(2);

	  shapes[0] = (1 - x) * (1 - y) * (1 - z);
	  shapes[1] = x * (1 - y) * (1 - z);
	  shapes[2] = x * y * (1 - z);
	  shapes[3] = (1 - x) * y * (1 - z);
	  shapes[4] = (1 - x) * (1 - y) * (z);
	  shapes[5] = x * (1 - y) * (z);
	  shapes[6] = x * y * (z);
	  shapes[7] = (1 - x) * y * (z);

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
			  T lame = shapes[e[i][0]] + shapes[e[i][1]];
			  T xi = sigma[e[i][1]] - sigma[e[i][0]];
			  shapes[8 + i] = (1 - xi * xi) * lame;
		  }
		  for (int i = 0; i < 12; i++)
		  {
			  shapes[e[i][0]] -= 0.5 * shapes[8 + i];
			  shapes[e[i][1]] -= 0.5 * shapes[8 + i];
		  }
		  break;
	  }

	  default:
	throw new Exception("CurvedElements::CalcShape 3d, element type not handled");

	};
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetCoefficients(ElementInfo & info, Vec<3> * coefs) const
  private void GetCoefficients(ElementInfo info, Vec < 3>[] coefs)
  {
	Element el = mesh[info.elnr];

	for (int i = 0; i < info.nv; i++)
	{
	  coefs[i] = Vec < 3> (mesh[el[i]]);
	}

	if (info.order == 1)
	{
		return;
	}

	int ii = info.nv;

	for (int i = 0; i < info.nedges; i++)
	{
	int first = edgecoeffsindex[info.edgenrs[i]];
	int next = edgecoeffsindex[info.edgenrs[i] + 1];
	for (int j = first; j < next; j++, ii++)
	{
	  coefs[ii] = edgecoeffs[j];
	}
	}
	for (int i = 0; i < info.nfaces; i++)
	{
	int first = facecoeffsindex[info.facenrs[i]];
	int next = facecoeffsindex[info.facenrs[i] + 1];
	for (int j = first; j < next; j++, ii++)
	{
	  coefs[ii] = facecoeffs[j];
	}
	}
  }

//C++ TO C# CONVERTER TODO TASK: The original C++ template specifier was replaced with a C# generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <typename T>
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void CalcElementDShapes(ElementInfo & info, const Point<3,T> xi, MatrixFixWidth<3,T> dshapes) const
  private void CalcElementDShapes<T>(ElementInfo info, Point<3,T> xi, MatrixFixWidth<3,T> dshapes)
  {
	// static int timer = NgProfiler::CreateTimer ("calcelementdshapes");

	Element el = mesh[info.elnr];

	// dshapes.SetSize(info.ndof);
//C++ TO C# CONVERTER TODO TASK: There is no equivalent in C# to 'alignof':
	if ((int)(dshapes(0,0)) % alignof(T) != 0)
	{
	  throw new Exception("alignment problem");
	}
	if (dshapes.Height() != info.ndof)
	{
	  throw new Exception("wrong height");
	}
	if (rational && info.order >= 2)
	{
	T w = 1;
	T[] dw = {0, 0, 0};

	T[] lami = {xi(0), xi(1), xi(2), 1 - xi(0) - xi(1) - xi(2)};
	T[][] dlami =
	{
		new T[] {1, 0, 0},
		new T[] {0, 1, 0},
		new T[] {0, 0, 1},
		new T[] {-1, -1, -1}
	};
	T[] shapes = Arrays.InitializeWithDefaultInstances<T>(10);

	for (int j = 0; j < 4; j++)
	{
		shapes[j] = lami[j] * lami[j];
		dshapes(j,0) = 2 * lami[j] * dlami[j][0];
		dshapes(j,1) = 2 * lami[j] * dlami[j][1];
		dshapes(j,2) = 2 * lami[j] * dlami[j][2];
	}

	ELEMENT_EDGE edges = MeshTopology.GetEdges1(TET);
	for (int j = 0; j < 6; j++)
	{
		T wi = edgeweight[info.edgenrs[j]];

		shapes[j + 4] = 2 * wi * lami[edges[j][0] - 1] * lami[edges[j][1] - 1];
		for (int k = 0; k < 3; k++)
		{
		  dshapes(j + 4,k) = 2 * wi * (lami[edges[j][0] - 1] * dlami[edges[j][1] - 1][k] + lami[edges[j][1] - 1] * dlami[edges[j][0] - 1][k]);
		}

		w += (wi - 1) * 2 * lami[edges[j][0] - 1] * lami[edges[j][1] - 1];
		for (int k = 0; k < 3; k++)
		{
		  dw[k] += 2 * (wi - 1) * (lami[edges[j][0] - 1] * dlami[edges[j][1] - 1][k] + lami[edges[j][1] - 1] * dlami[edges[j][0] - 1][k]);
		}
	}
	// shapes *= 1.0 / w;
	dshapes *= 1.0 / w;
	for (int i = 0; i < 10; i++)
	{
	  for (int j = 0; j < 3; j++)
	  {
		dshapes(i,j) -= shapes[i] * dw[j] / (w * w);
	  }
	}
	return;
	}

	/*
	if (typeid(T) == typeid(SIMD<double>))
	  {
	    if (el.GetType() == HEX)
	      dshapes = T(0.0);
	    return;
	  }
	*/
	switch (el.GetType())
	{
	  case TET:
	  {
		  // if (typeid(T) == typeid(SIMD<double>)) return;

		  dshapes = T(0.0);

	  dshapes(0,0) = 1;
	  dshapes(1,1) = 1;
	  dshapes(2,2) = 1;
	  dshapes(3,0) = -1;
	  dshapes(3,1) = -1;
	  dshapes(3,2) = -1;

	  if (info.order == 1)
	  {
		  return;
	  }

	  T[] lami = {xi(0), xi(1), xi(2), 1 - xi(0) - xi(1) - xi(2)};
	  int ii = 4;
	  ELEMENT_EDGE edges = MeshTopology.GetEdges1(TET);
	  for (int i = 0; i < 6; i++)
	  {
		  int eorder = edgeorder[info.edgenrs[i]];
		  if (eorder >= 2)
		  {
		  int vi1 = edges[i][0] - 1;
		  int vi2 = edges[i][1] - 1;
		  if (el[vi1] > el[vi2])
		  {
			  swap(vi1, vi2);
		  }

		  netgen.GlobalMembers.CalcScaledEdgeShapeDxDt < 3> (eorder, lami[vi1] - lami[vi2], lami[vi1] + lami[vi2], dshapes(ii,0));

		  Mat<2,3,T> trans = new Mat<2,3,T>();
		  for (int j = 0; j < 3; j++)
		  {
			  trans(0,j) = dshapes(vi1,j) - dshapes(vi2,j);
			  trans(1,j) = dshapes(vi1,j) + dshapes(vi2,j);
		  }

		  for (int j = 0; j < order - 1; j++)
		  {
			  T ddx = dshapes(ii + j,0);
			  T ddt = dshapes(ii + j,1);
			  dshapes(ii + j,0) = ddx * trans(0,0) + ddt * trans(1,0);
			  dshapes(ii + j,1) = ddx * trans(0,1) + ddt * trans(1,1);
			  dshapes(ii + j,2) = ddx * trans(0,2) + ddt * trans(1,2);
		  }

		  ii += eorder - 1;
		  }
	  }

	  ELEMENT_FACE faces = MeshTopology.GetFaces1(TET);
	  for (int i = 0; i < 4; i++)
	  {
		  int forder = faceorder[info.facenrs[i]];
		  if (forder >= 3)
		  {
		  int[] fnums = {faces[i][0] - 1, faces[i][1] - 1, faces[i][2] - 1};
		  if (el[fnums[0]] > el[fnums[1]])
		  {
			  swap(fnums[0], fnums[1]);
		  }
		  if (el[fnums[1]] > el[fnums[2]])
		  {
			  swap(fnums[1], fnums[2]);
		  }
		  if (el[fnums[0]] > el[fnums[1]])
		  {
			  swap(fnums[0], fnums[1]);
		  }

		  netgen.GlobalMembers.CalcScaledTrigShapeDxDyDt(forder, lami[fnums[1]] - lami[fnums[0]], lami[fnums[2]], lami[fnums[0]] + lami[fnums[1]] + lami[fnums[2]], dshapes(ii,0));

		  Mat<3,3,T> trans = new Mat<3,3,T>();
		  for (int j = 0; j < 3; j++)
		  {
			  trans(0,j) = dshapes(fnums[1],j) - dshapes(fnums[0],j);
			  trans(1,j) = dshapes(fnums[2],j);
			  trans(2,j) = dshapes(fnums[0],j) + dshapes(fnums[1],j) + dshapes(fnums[2],j);
		  }

		  int nfd = (forder - 1) * (forder - 2) / 2;
		  for (int j = 0; j < nfd; j++)
		  {
			  T ddx = dshapes(ii + j,0);
			  T ddy = dshapes(ii + j,1);
			  T ddt = dshapes(ii + j,2);
			  dshapes(ii + j,0) = ddx * trans(0,0) + ddy * trans(1,0) + ddt * trans(2,0);
			  dshapes(ii + j,1) = ddx * trans(0,1) + ddy * trans(1,1) + ddt * trans(2,1);
			  dshapes(ii + j,2) = ddx * trans(0,2) + ddy * trans(1,2) + ddt * trans(2,2);
		  }

		  ii += nfd;
		  }
	  }

	  break;
	  }

	  case TET10:
	  {
		  // if (typeid(T) == typeid(SIMD<double>)) return;

	  if (dshapes.Height() == 4)
	  {
		  dshapes = T(0.0);

		  dshapes(0,0) = 1;
		  dshapes(1,1) = 1;
		  dshapes(2,2) = 1;
		  dshapes(3,0) = -1;
		  dshapes(3,1) = -1;
		  dshapes(3,2) = -1;
	  }
	  else
	  {
		  AutoDiff<3,T> x = new AutoDiff<3,T>(xi(0), 0);
		  AutoDiff<3,T> y = new AutoDiff<3,T>(xi(1), 1);
		  AutoDiff<3,T> z = new AutoDiff<3,T>(xi(2), 2);
		  AutoDiff<3,T> lam4 = 1 - x - y - z;
		  AutoDiff<3,T>[] shapes = Arrays.InitializeWithDefaultInstances<AutoDiff>(10);

		  shapes[0] = 2 * x * x - x;
		  shapes[1] = 2 * y * y - y;
		  shapes[2] = 2 * z * z - z;
		  shapes[3] = 2 * lam4 * lam4 - lam4;

		  shapes[4] = 4 * x * y;
		  shapes[5] = 4 * x * z;
		  shapes[6] = 4 * x * lam4;
		  shapes[7] = 4 * y * z;
		  shapes[8] = 4 * y * lam4;
		  shapes[9] = 4 * z * lam4;

		  for (int i = 0; i < 10; i++)
		  {
		  dshapes(i,0) = shapes[i].DValue(0);
		  dshapes(i,1) = shapes[i].DValue(1);
		  dshapes(i,2) = shapes[i].DValue(2);
		  }

	  }
	  break;

	  break;
	  }


	  case PRISM:
	  {
	  T[] lami = {xi(0), xi(1), 1 - xi(0) - xi(1), xi(0), xi(1), 1 - xi(0) - xi(1)};
	  T[] lamiz = {1 - xi(2), 1 - xi(2), 1 - xi(2), xi(2), xi(2), xi(2)};
	  T[] dlamiz = {-1, -1, -1, 1, 1, 1};
	  T[][] dlami =
	  {
		  new T[] {1, 0},
		  new T[] {0, 1},
		  new T[] {-1, -1},
		  new T[] {1, 0},
		  new T[] {0, 1},
		  new T[] {-1, -1}
	  };
	  for (int i = 0; i < 6; i++)
	  {
		  // shapes(i) = lami[i%3] * ( (i < 3) ? (1-xi(2)) : xi(2) );
		  dshapes(i,0) = dlami[i % 3][0] * ((i < 3) ? (1 - xi(2)) : xi(2));
		  dshapes(i,1) = dlami[i % 3][1] * ((i < 3) ? (1 - xi(2)) : xi(2));
		  dshapes(i,2) = lami[i % 3] * ((i < 3) ? -1 : 1);
	  }

	  int ii = 6;

	  if (info.order == 1)
	  {
		  return;
	  }


	  ELEMENT_EDGE edges = MeshTopology.GetEdges1(PRISM);
	  for (int i = 0; i < 6; i++) // horizontal edges
	  {
		  int order = edgeorder[info.edgenrs[i]];
		  if (order >= 2)
		  {
		  int vi1 = (edges[i][0] - 1);
		  int vi2 = (edges[i][1] - 1);
		  if (el[vi1] > el[vi2])
		  {
			  swap(vi1, vi2);
		  }
		  vi1 = vi1 % 3;
		  vi2 = vi2 % 3;

				  ArrayMem<T,20> shapei_mem = new ArrayMem<T,20>((uint)(order + 1));
		  TFlatVector<T> shapei = new TFlatVector<T>(order + 1, shapei_mem[0]);
		  netgen.GlobalMembers.CalcScaledEdgeShapeDxDt < 3> (order, lami[vi1] - lami[vi2], lami[vi1] + lami[vi2], dshapes(ii,0));
		  netgen.GlobalMembers.CalcScaledEdgeShape(order, lami[vi1] - lami[vi2], lami[vi1] + lami[vi2], shapei(0));

		  Mat<2,2,T> trans = new Mat<2,2,T>();
		  for (int j = 0; j < 2; j++)
		  {
			  trans(0,j) = dlami[vi1][j] - dlami[vi2][j];
			  trans(1,j) = dlami[vi1][j] + dlami[vi2][j];
		  }

		  for (int j = 0; j < order - 1; j++)
		  {
			  T ddx = dshapes(ii + j,0);
			  T ddt = dshapes(ii + j,1);
			  dshapes(ii + j,0) = ddx * trans(0,0) + ddt * trans(1,0);
			  dshapes(ii + j,1) = ddx * trans(0,1) + ddt * trans(1,1);
		  }



		  T facz = (i < 3) ? (1 - xi(2)) : xi(2);
		  T dfacz = (i < 3) ? (-1) : 1;
		  for (int j = 0; j < order - 1; j++)
		  {
			  dshapes(ii + j,0) *= facz;
			  dshapes(ii + j,1) *= facz;
			  dshapes(ii + j,2) = shapei(j) * dfacz;
		  }

		  ii += order - 1;
		  }
	  }

		  // if (typeid(T) == typeid(SIMD<double>)) return;


	  for (int i = 6; i < 9; i++) // vertical edges
	  {
		  int eorder = edgeorder[info.edgenrs[i]];
		  if (eorder >= 2)
		  {
		  int vi1 = (edges[i][0] - 1);
		  int vi2 = (edges[i][1] - 1);
		  if (el[vi1] > el[vi2])
		  {
			  swap(vi1, vi2);
		  }

		  T bubz = lamiz[vi1] * lamiz[vi2];
		  T dbubz = dlamiz[vi1] * lamiz[vi2] + lamiz[vi1] * dlamiz[vi2];
		  T polyz = lamiz[vi1] - lamiz[vi2];
		  T dpolyz = dlamiz[vi1] - dlamiz[vi2];
		  T bubxy = lami[(vi1) % 3];
		  T dbubxydx = dlami[(vi1) % 3][0];
		  T dbubxydy = dlami[(vi1) % 3][1];

		  for (int j = 0; j < eorder - 1; j++)
		  {
			  dshapes(ii + j,0) = dbubxydx * bubz;
			  dshapes(ii + j,1) = dbubxydy * bubz;
			  dshapes(ii + j,2) = bubxy * dbubz;

			  dbubz = bubz * dpolyz + dbubz * polyz;
			  bubz *= polyz;
		  }
		  ii += eorder - 1;
		  }
	  }


	  if (info.order == 2)
	  {
		  return;
	  }
	  // FACE SHAPES
	  ELEMENT_FACE faces = MeshTopology.GetFaces1(PRISM);
	  for (int i = 0; i < 2; i++)
	  {
		  int forder = faceorder[info.facenrs[i]];

		  if (forder < 3)
		  {
			  continue;
		  }
		  int ndf = (forder + 1) * (forder + 2) / 2 - 3 - 3 * (forder - 1);

		  int[] fav = {faces[i][0] - 1, faces[i][1] - 1, faces[i][2] - 1};
		  if (el[fav[0]] > el[fav[1]])
		  {
			  swap(fav[0],fav[1]);
		  }
		  if (el[fav[1]] > el[fav[2]])
		  {
			  swap(fav[1],fav[2]);
		  }
		  if (el[fav[0]] > el[fav[1]])
		  {
			  swap(fav[0],fav[1]);
		  }

			  ArrayMem<T,2 20> dshapei_mem = new ArrayMem<T,2 20>((uint)ndf);
			  ArrayMem<T,20> shapei_mem = new ArrayMem<T,20>((uint)ndf);
		  MatrixFixWidth<2,T> dshapei = new MatrixFixWidth<2,T>(ndf, dshapei_mem[0]);
		  TFlatVector<T> shapei = new TFlatVector<T>(ndf, shapei_mem[0]);

		  netgen.GlobalMembers.CalcTrigShapeDxDy(forder, lami[fav[2]] - lami[fav[1]], lami[fav[0]], dshapei(0,0));
		  netgen.GlobalMembers.CalcTrigShape(forder, lami[fav[2]] - lami[fav[1]], lami[fav[0]], shapei(0));

		  Mat<2,2,T> trans = new Mat<2,2,T>();
		  for (int j = 0; j < 2; j++)
		  {
		  trans(0,j) = dlami[fav[2]][j] - dlami[fav[1]][j];
		  trans(1,j) = dlami[fav[0]][j];
		  }

		  for (int j = 0; j < ndf; j++)
		  {
		  // double ddx = dshapes(ii+j,0);
		  // double ddt = dshapes(ii+j,1);
		  T ddx = dshapei(j,0);
		  T ddt = dshapei(j,1);
		  dshapes(ii + j,0) = ddx * trans(0,0) + ddt * trans(1,0);
		  dshapes(ii + j,1) = ddx * trans(0,1) + ddt * trans(1,1);
		  }

		  for (int j = 0; j < ndf; j++)
		  {
		  dshapes(ii + j,0) *= lamiz[fav[1]];
		  dshapes(ii + j,1) *= lamiz[fav[1]];
		  dshapes(ii + j,2) = shapei(j) * dlamiz[fav[1]];
		  }
		  ii += ndf;
	  }

	  break;

	  }

	  case PRISM15:
	  {
		  AutoDiff<3,T> x = new AutoDiff<3,T>(xi(0), 0);
		  AutoDiff<3,T> y = new AutoDiff<3,T>(xi(1), 1);
		  AutoDiff<3,T> z = new AutoDiff<3,T>(xi(2), 2);
		  AutoDiff<3,T>[] ad = Arrays.InitializeWithDefaultInstances<AutoDiff>(15);
		  AutoDiff<3,T> lam = 1 - x - y;
		  AutoDiff<3,T> lamz = 1 - z;

		  ad[0] = (2 * x * x - x) * (2 * lamz * lamz - lamz);
		  ad[1] = (2 * y * y - y) * (2 * lamz * lamz - lamz);
		  ad[2] = (2 * lam * lam - lam) * (2 * lamz * lamz - lamz);
		  ad[3] = (2 * x * x - x) * (2 * z * z - z);
		  ad[4] = (2 * y * y - y) * (2 * z * z - z);
		  ad[5] = (2 * lam * lam - lam) * (2 * z * z - z);
		  ad[6] = 4 * x * y * (2 * lamz * lamz - lamz);
		  ad[7] = 4 * x * lam * (2 * lamz * lamz - lamz);
		  ad[8] = 4 * y * lam * (2 * lamz * lamz - lamz);
		  ad[9] = x * 4 * z * (1 - z);
		  ad[10] = y * 4 * z * (1 - z);
		  ad[11] = lam * 4 * z * (1 - z);
		  ad[12] = 4 * x * y * (2 * z * z - z);
		  ad[13] = 4 * x * lam * (2 * z * z - z);
		  ad[14] = 4 * y * lam * (2 * z * z - z);

		  for (int i = 0; i < 15; i++)
		  {
			for (int j = 0; j < 3; j++)
			{
			  dshapes(i,j) = ad[i].DValue(j);
			}
		  }
		  break;
	  }
	  case PYRAMID:
	  {
		  // if (typeid(T) == typeid(SIMD<double>)) return;

	  dshapes = T(0.0);
	  T x = xi(0);
	  T y = xi(1);
	  T z = xi(2);

	  // if (z == 1.) z = 1-1e-10;
		  z *= 1 - 1e-12;
	  T z1 = 1 - z;
	  T z2 = z1 * z1;

	  dshapes(0,0) = -(z1 - y) / z1;
	  dshapes(0,1) = -(z1 - x) / z1;
	  dshapes(0,2) = ((x + y + 2 * z - 2) * z1 + (z1 - y) * (z1 - x)) / z2;

	  dshapes(1,0) = (z1 - y) / z1;
	  dshapes(1,1) = -x / z1;
	  dshapes(1,2) = (-x * z1 + x * (z1 - y)) / z2;

	  dshapes(2,0) = y / z1;
	  dshapes(2,1) = x / z1;
	  dshapes(2,2) = x * y / z2;

	  dshapes(3,0) = -y / z1;
	  dshapes(3,1) = (z1 - x) / z1;
	  dshapes(3,2) = (-y * z1 + y * (z1 - x)) / z2;

	  dshapes(4,0) = 0;
	  dshapes(4,1) = 0;
	  dshapes(4,2) = 1;

	  if (info.order == 1)
	  {
		  return;
	  }

	  int ii = 5;
	  ELEMENT_EDGE edges = MeshTopology.GetEdges1(PYRAMID);
	  // if (z == 1.) z = 1-1e-10;
		  z *= 1 - 1e-12;
		  T[] shapes = Arrays.InitializeWithDefaultInstances<T>(5);
	  shapes[0] = (1 - z - x) * (1 - z - y) / (1 - z);
	  shapes[1] = x * (1 - z - y) / (1 - z);
	  shapes[2] = x * y / (1 - z);
	  shapes[3] = (1 - z - x) * y / (1 - z);
	  shapes[4] = z;

		  T[] sigma = {((1 - z - x) + (1 - z - y)), (x + (1 - z - y)), (x + y), ((1 - z - x) + y)};
		  T[][] dsigma =
		  {
			  new T[] {-1, -1, -2},
			  new T[] {1, -1, -1},
			  new T[] {1, 1, 0},
			  new T[] {-1, 1, -1}
		  };
		  T[] dz = {0, 0, 1};
	  for (int i = 0; i < 4; i++) // horizontal edges
	  {
		  int eorder = edgeorder[info.edgenrs[i]];
		  if (eorder >= 2)
		  {
		  int vi1 = (edges[i][0] - 1);
		  int vi2 = (edges[i][1] - 1);
		  if (el[vi1] > el[vi2])
		  {
			  swap(vi1, vi2);
		  }

				  ArrayMem<T,20> shapei_mem = new ArrayMem<T,20>((uint)(eorder + 1));
		  TFlatVector<T> shapei = new TFlatVector<T>(eorder + 1, shapei_mem[0]);
		  netgen.GlobalMembers.CalcScaledEdgeShapeDxDt < 3> (eorder, sigma[vi1] - sigma[vi2], 1 - z, dshapes(ii,0));
		  netgen.GlobalMembers.CalcScaledEdgeShape(eorder, sigma[vi1] - sigma[vi2], 1 - z, shapei(0));
		  T fac = (shapes[vi1] + shapes[vi2]) / (1 - z);
				  T[] dfac = Arrays.InitializeWithDefaultInstances<T>(3);
				  for (int k = 0; k < 3; k++)
				  {
					dfac[k] = ((dshapes(vi1,k) + dshapes(vi2,k)) * (1 - z) - (shapes[vi1] + shapes[vi2]) * (-dshapes(4,k))) / netgen.GlobalMembers.sqr(1 - z);
				  }

		  for (int j = 0; j < eorder - 1; j++)
		  {
					  T ddx = dshapes(ii + j,0);
					  T ddt = dshapes(ii + j,1);
					  for (int k = 0; k < 3; k++)
					  {
						dshapes(ii + j,k) = fac * (ddx * (dsigma[vi1][k] - dsigma[vi2][k]) - ddt * dz[k]) + dfac[k] * shapei(j);
					  }
		  }

		  ii += eorder - 1;
		  }
	  }

	  break;
	  }

	  case PYRAMID13:
	  {
	  T x = xi(0);
	  T y = xi(1);
	  T z = xi(2);
		  z *= 1 - 1e-12;
		  dshapes(0,0) = 0.5 * z - 0.5 * z * (2 * x + z - 1) * (2 * y + z - 1) / (-z + 1) - 0.5 * (-2 * x - z + 2) * (-2 * y - z + 2) + (-0.5 * x - 0.5 * y - 0.5 * z + 0.25) * (4 * y + 2 * z + 2 * z * (2 * y + z - 1) / (-z + 1) - 4);
		  dshapes(0,1) = 0.5 * z - 0.5 * z * (2 * x + z - 1) * (2 * y + z - 1) / (-z + 1) - 0.5 * (-2 * x - z + 2) * (-2 * y - z + 2) + (-0.5 * x - 0.5 * y - 0.5 * z + 0.25) * (4 * x + 2 * z + 2 * z * (2 * x + z - 1) / (-z + 1) - 4);
		  dshapes(0,2) = 0.5 * z - 0.5 * z * (2 * x + z - 1) * (2 * y + z - 1) / (-z + 1) - 0.5 * (-2 * x - z + 2) * (-2 * y - z + 2) + (-0.5 * x - 0.5 * y - 0.5 * z + 0.25) * (2 * x + 2 * y + 2 * z + z * (2 * x + z - 1) / (-z + 1) + z * (2 * y + z - 1) / (-z + 1) + z * (2 * x + z - 1) * (2 * y + z - 1) / ((-z + 1) * (-z + 1)) - 5 + (2 * x + z - 1) * (2 * y + z - 1) / (-z + 1));
		  dshapes(1,0) = -0.5 * z - 0.5 * z * (2 * x + z - 1) * (2 * y + z - 1) / (-z + 1) + 0.5 * (2 * x + z) * (-2 * y - z + 2) + (0.5 * x - 0.5 * y - 0.25) * (-4 * y - 2 * z - 2 * z * (2 * y + z - 1) / (-z + 1) + 4);
		  dshapes(1,1) = 0.5 * z + 0.5 * z * (2 * x + z - 1) * (2 * y + z - 1) / (-z + 1) - 0.5 * (2 * x + z) * (-2 * y - z + 2) + (-4 * x - 2 * z - 2 * z * (2 * x + z - 1) / (-z + 1)) * (0.5 * x - 0.5 * y - 0.25);
		  dshapes(1,2) = (0.5 * x - 0.5 * y - 0.25) * (-2 * x - 2 * y - 2 * z - z * (2 * x + z - 1) / (-z + 1) - z * (2 * y + z - 1) / (-z + 1) - z * (2 * x + z - 1) * (2 * y + z - 1) / ((-z + 1) * (-z + 1)) + 1 - (2 * x + z - 1) * (2 * y + z - 1) / (-z + 1));
		  dshapes(2,0) = -0.5 * z + 0.5 * z * (2 * x + z - 1) * (2 * y + z - 1) / (-z + 1) + 0.5 * (2 * x + z) * (2 * y + z) + (4 * y + 2 * z + 2 * z * (2 * y + z - 1) / (-z + 1)) * (0.5 * x + 0.5 * y + 0.5 * z - 0.75);
		  dshapes(2,1) = -0.5 * z + 0.5 * z * (2 * x + z - 1) * (2 * y + z - 1) / (-z + 1) + 0.5 * (2 * x + z) * (2 * y + z) + (4 * x + 2 * z + 2 * z * (2 * x + z - 1) / (-z + 1)) * (0.5 * x + 0.5 * y + 0.5 * z - 0.75);
		  dshapes(2,2) = -0.5 * z + 0.5 * z * (2 * x + z - 1) * (2 * y + z - 1) / (-z + 1) + 0.5 * (2 * x + z) * (2 * y + z) + (0.5 * x + 0.5 * y + 0.5 * z - 0.75) * (2 * x + 2 * y + 2 * z + z * (2 * x + z - 1) / (-z + 1) + z * (2 * y + z - 1) / (-z + 1) + z * (2 * x + z - 1) * (2 * y + z - 1) / ((-z + 1) * (-z + 1)) - 1 + (2 * x + z - 1) * (2 * y + z - 1) / (-z + 1));
		  dshapes(3,0) = 0.5 * z + 0.5 * z * (2 * x + z - 1) * (2 * y + z - 1) / (-z + 1) - 0.5 * (2 * y + z) * (-2 * x - z + 2) + (-0.5 * x + 0.5 * y - 0.25) * (-4 * y - 2 * z - 2 * z * (2 * y + z - 1) / (-z + 1));
		  dshapes(3,1) = -0.5 * z - 0.5 * z * (2 * x + z - 1) * (2 * y + z - 1) / (-z + 1) + 0.5 * (2 * y + z) * (-2 * x - z + 2) + (-0.5 * x + 0.5 * y - 0.25) * (-4 * x - 2 * z - 2 * z * (2 * x + z - 1) / (-z + 1) + 4);
		  dshapes(3,2) = (-0.5 * x + 0.5 * y - 0.25) * (-2 * x - 2 * y - 2 * z - z * (2 * x + z - 1) / (-z + 1) - z * (2 * y + z - 1) / (-z + 1) - z * (2 * x + z - 1) * (2 * y + z - 1) / ((-z + 1) * (-z + 1)) + 1 - (2 * x + z - 1) * (2 * y + z - 1) / (-z + 1));
		  dshapes(4,0) = 0;
		  dshapes(4,1) = 0;
		  dshapes(4,2) = 4 * z - 1;
		  dshapes(5,0) = -4 * x * (-2 * y - 2 * z + 2) / (-2 * z + 2) + 2 * (-2 * x - 2 * z + 2) * (-2 * y - 2 * z + 2) / (-2 * z + 2);
		  dshapes(5,1) = -4 * x * (-2 * x - 2 * z + 2) / (-2 * z + 2);
		  dshapes(5,2) = -4 * x * (-2 * x - 2 * z + 2) / (-2 * z + 2) - 4 * x * (-2 * y - 2 * z + 2) / (-2 * z + 2) + 4 * x * (-2 * x - 2 * z + 2) * (-2 * y - 2 * z + 2) / ((-2 * z + 2) * (-2 * z + 2));
		  dshapes(6,0) = -8 * x * y / (-2 * z + 2) + 4 * y * (-2 * x - 2 * z + 2) / (-2 * z + 2);
		  dshapes(6,1) = 4 * x * (-2 * x - 2 * z + 2) / (-2 * z + 2);
		  dshapes(6,2) = -8 * x * y / (-2 * z + 2) + 8 * x * y * (-2 * x - 2 * z + 2) / ((-2 * z + 2) * (-2 * z + 2));
		  dshapes(7,0) = -4 * y * (-2 * y - 2 * z + 2) / (-2 * z + 2);
		  dshapes(7,1) = -4 * y * (-2 * x - 2 * z + 2) / (-2 * z + 2) + 2 * (-2 * x - 2 * z + 2) * (-2 * y - 2 * z + 2) / (-2 * z + 2);
		  dshapes(7,2) = -4 * y * (-2 * x - 2 * z + 2) / (-2 * z + 2) - 4 * y * (-2 * y - 2 * z + 2) / (-2 * z + 2) + 4 * y * (-2 * x - 2 * z + 2) * (-2 * y - 2 * z + 2) / ((-2 * z + 2) * (-2 * z + 2));
		  dshapes(8,0) = 4 * y * (-2 * y - 2 * z + 2) / (-2 * z + 2);
		  dshapes(8,1) = -8 * x * y / (-2 * z + 2) + 4 * x * (-2 * y - 2 * z + 2) / (-2 * z + 2);
		  dshapes(8,2) = -8 * x * y / (-2 * z + 2) + 8 * x * y * (-2 * y - 2 * z + 2) / ((-2 * z + 2) * (-2 * z + 2));
		  dshapes(9,0) = -2 * z * (-2 * y - 2 * z + 2) / (-z + 1);
		  dshapes(9,1) = -2 * z * (-2 * x - 2 * z + 2) / (-z + 1);
		  dshapes(9,2) = -2 * z * (-2 * x - 2 * z + 2) / (-z + 1) - 2 * z * (-2 * y - 2 * z + 2) / (-z + 1) + z * (-2 * x - 2 * z + 2) * (-2 * y - 2 * z + 2) / ((-z + 1) * (-z + 1)) + (-2 * x - 2 * z + 2) * (-2 * y - 2 * z + 2) / (-z + 1);
		  dshapes(10,0) = 2 * z * (-2 * y - 2 * z + 2) / (-z + 1);
		  dshapes(10,1) = -4 * x * z / (-z + 1);
		  dshapes(10,2) = -4 * x * z / (-z + 1) + 2 * x * z * (-2 * y - 2 * z + 2) / ((-z + 1) * (-z + 1)) + 2 * x * (-2 * y - 2 * z + 2) / (-z + 1);
		  dshapes(11,0) = 4 * y * z / (-z + 1);
		  dshapes(11,1) = 4 * x * z / (-z + 1);
		  dshapes(11,2) = 4 * x * y * z / ((-z + 1) * (-z + 1)) + 4 * x * y / (-z + 1);
		  dshapes(12,0) = -4 * y * z / (-z + 1);
		  dshapes(12,1) = 2 * z * (-2 * x - 2 * z + 2) / (-z + 1);
		  dshapes(12,2) = -4 * y * z / (-z + 1) + 2 * y * z * (-2 * x - 2 * z + 2) / ((-z + 1) * (-z + 1)) + 2 * y * (-2 * x - 2 * z + 2) / (-z + 1);
		  break;
	  }

	  case HEX:
	  {
		  // if (typeid(T) == typeid(SIMD<double>)) return;

		  // NgProfiler::StartTimer(timer);
	  T x = xi(0);
	  T y = xi(1);
	  T z = xi(2);

	  // shapes[0] = (1-x)*(1-y)*(1-z);
	  dshapes(0,0) = - (1 - y) * (1 - z);
	  dshapes(0,1) = (1 - x) * (-1) * (1 - z);
	  dshapes(0,2) = (1 - x) * (1 - y) * (-1);

	  // shapes[1] =    x *(1-y)*(1-z);
	  dshapes(1,0) = (1 - y) * (1 - z);
	  dshapes(1,1) = -x * (1 - z);
	  dshapes(1,2) = -x * (1 - y);

	  // shapes[2] =    x *   y *(1-z);
	  dshapes(2,0) = y * (1 - z);
	  dshapes(2,1) = x * (1 - z);
	  dshapes(2,2) = -x * y;

	  // shapes[3] = (1-x)*   y *(1-z);
	  dshapes(3,0) = -y * (1 - z);
	  dshapes(3,1) = (1 - x) * (1 - z);
	  dshapes(3,2) = -(1 - x) * y;

	  // shapes[4] = (1-x)*(1-y)*z;
	  dshapes(4,0) = - (1 - y) * z;
	  dshapes(4,1) = (1 - x) * (-1) * z;
	  dshapes(4,2) = (1 - x) * (1 - y) * 1;

	  // shapes[5] =    x *(1-y)*z;
	  dshapes(5,0) = (1 - y) * z;
	  dshapes(5,1) = -x * z;
	  dshapes(5,2) = x * (1 - y);

	  // shapes[6] =    x *   y *z;
	  dshapes(6,0) = y * z;
	  dshapes(6,1) = x * z;
	  dshapes(6,2) = x * y;

	  // shapes[7] = (1-x)*   y *z;
	  dshapes(7,0) = -y * z;
	  dshapes(7,1) = (1 - x) * z;
	  dshapes(7,2) = (1 - x) * y;

		  // NgProfiler::StopTimer(timer);

	  if (info.order == 1)
	  {
		  return;
	  }

	  T[] shapes = {(1 - x) * (1 - y) * (1 - z), x * (1 - y) * (1 - z), x * y * (1 - z), (1 - x) * y * (1 - z), (1 - x) * (1 - y) * (z), x * (1 - y) * (z), x * y * (z), (1 - x) * y * (z)};

	  T[] mu = {(1 - x) + (1 - y) + (1 - z), x + (1 - y) + (1 - z), x + y + (1 - z), (1 - x) + y + (1 - z), (1 - x) + (1 - y) + (z), x + (1 - y) + (z), x + y + (z), (1 - x) + y + (z)};

	  T[][] dmu =
	  {
		  new T[] {-1, -1, -1},
		  new T[] {1, -1, -1},
		  new T[] {1, 1, -1},
		  new T[] {-1, 1, -1},
		  new T[] {-1, -1, 1},
		  new T[] {1, -1, 1},
		  new T[] {1, 1, 1},
		  new T[] {-1, 1, 1}
	  };

	  ArrayMem<T, 20> hshapes = new ArrayMem<T, 20>(order + 1);
	  ArrayMem<T, 20> hdshapes = new ArrayMem<T, 20>(order + 1);

	  int ii = 8;
	  ELEMENT_EDGE edges = MeshTopology.GetEdges1(HEX);

	  for (int i = 0; i < 8; i++)
	  {
		  int eorder = edgeorder[info.edgenrs[i]];
		  if (eorder >= 2)
		  {
		  int vi1 = edges[i][0] - 1;
		  int vi2 = edges[i][1] - 1;
		  if (el[vi1] > el[vi2])
		  {
			  swap(vi1, vi2);
		  }

		  netgen.GlobalMembers.CalcEdgeShapeDx(eorder, mu[vi1] - mu[vi2], hshapes[0], hdshapes[0]);

		  T lame = shapes[vi1] + shapes[vi2];
		  T[] dlame = {dshapes(vi1, 0) + dshapes(vi2, 0), dshapes(vi1, 1) + dshapes(vi2, 1), dshapes(vi1, 2) + dshapes(vi2, 2)};

		  for (int j = 0; j < eorder - 1; j++)
		  {
			for (int k = 0; k < 3; k++)
			{
			  dshapes(ii + j, k) = lame * hdshapes[j] * (dmu[vi1][k] - dmu[vi2][k]) + dlame[k] * hshapes[j];
			}
		  }

		  ii += eorder - 1;
		  }
	  }

	  /*
	   *testout << "quad, dshape = " << endl << dshapes << endl;
	   for (int i = 0; i < 2; i++)
	   {
	   Point<2> xil = xi, xir = xi;
	   Vector shapesl(dshapes.Height()), shapesr(dshapes.Height());
	   xil(i) -= 1e-6;
	   xir(i) += 1e-6;
	   CalcElementShapes (info, xil, shapesl);
	   CalcElementShapes (info, xir, shapesr);
	  
	   for (int j = 0; j < dshapes.Height(); j++)
	   dshapes(j,i) = 1.0 / 2e-6 * (shapesr(j)-shapesl(j));
	   }
	  
	   *testout << "quad, num dshape = " << endl << dshapes << endl;
	   */
	  break;
	  }
	  case HEX20:
	  {
		  AutoDiff<3,T> x = new AutoDiff<3,T>(xi(0), 0);
		  AutoDiff<3,T> y = new AutoDiff<3,T>(xi(1), 1);
		  AutoDiff<3,T> z = new AutoDiff<3,T>(xi(2), 2);
		  AutoDiff<3,T>[] ad = Arrays.InitializeWithDefaultInstances<AutoDiff>(20);

		  ad[0] = (1 - x) * (1 - y) * (1 - z);
	  ad[1] = x * (1 - y) * (1 - z);
	  ad[2] = x * y * (1 - z);
	  ad[3] = (1 - x) * y * (1 - z);
	  ad[4] = (1 - x) * (1 - y) * (z);
	  ad[5] = x * (1 - y) * (z);
	  ad[6] = x * y * (z);
	  ad[7] = (1 - x) * y * (z);

//C++ TO C# CONVERTER TODO TASK: The following line could not be converted:
		  AutoDiff<3,T> sigma[8] = {(1 - x) + (1 - y) + (1 - z),x + (1 - y) + (1 - z),x + y + (1 - z),(1 - x) + y + (1 - z), (1 - x) + (1 - y) + z,x + (1 - y) + z,x + y + z,(1 - x) + y + z};

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
			  var lame = ad[e[i][0]] + ad[e[i][1]];
			  var xi = sigma[e[i][1]] - sigma[e[i][0]];
			  ad[8 + i] = (1 - xi * xi) * lame;
		  }
		  for (int i = 0; i < 12; i++)
		  {
			  ad[e[i][0]] -= 0.5 * ad[8 + i];
			  ad[e[i][1]] -= 0.5 * ad[8 + i];
		  }
		  for (int i = 0; i < 20; i++)
		  {
			for (int j = 0; j < 3; j++)
			{
			  dshapes(i,j) = ad[i].DValue(j);
			}
		  }
		  break;
	  }
	  default:
	throw new Exception("CurvedElements::CalcDShape 3d, element type not handled");
	}

	/*
	  DenseMatrix dshapes2 (info.ndof, 3);
	  Vector shapesl(info.ndof);
	  Vector shapesr(info.ndof);
	
	  double eps = 1e-6;
	  for (int i = 0; i < 3; i++)
	  {
	  Point<3> xl = xi;
	  Point<3> xr = xi;

	  xl(i) -= eps;
	  xr(i) += eps;
	  CalcElementShapes (info, xl, shapesl);
	  CalcElementShapes (info, xr, shapesr);

	  for (int j = 0; j < info.ndof; j++)
	  dshapes2(j,i) = (shapesr(j)-shapesl(j)) / (2*eps);
	  }
	  (*testout) << "dshapes = " << endl << dshapes << endl;
	  (*testout) << "dshapes2 = " << endl << dshapes2 << endl;
	  dshapes2 -= dshapes;
	  (*testout) << "diff = " << endl << dshapes2 << endl;
	*/
  }

//C++ TO C# CONVERTER TODO TASK: The original C++ template specifier was replaced with a C# generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <typename T>

  // extern int mappingvar;
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool EvaluateMapping(ElementInfo & info, Point<3,T> xi, Point<3,T> & mx, Mat<3,3,T> & jac) const
  private bool EvaluateMapping<T>(ElementInfo info, Point<3,T> xi, Point<3,T> mx, Mat<3,3,T> jac)
  {
	Element el = mesh[info.elnr];
	if (rational && info.order >= 2)
	{
		return false; // not supported
	}

	AutoDiff<3,T> x = new AutoDiff<3,T>(xi(0), 0);
	AutoDiff<3,T> y = new AutoDiff<3,T>(xi(1), 1);
	AutoDiff<3,T> z = new AutoDiff<3,T>(xi(2), 2);

//C++ TO C# CONVERTER TODO TASK: The following line could not be converted:
	AutoDiff<3,T> mapped_x[3] = {T(0.0), T(0.0), T(0.0)};

	switch (el.GetType())
	{
	  case TET:
	  {
		  // if (info.order >= 2) return false; // not yet supported
//C++ TO C# CONVERTER TODO TASK: The following line could not be converted:
		  AutoDiff<3,T> lami[4] = {x, y, z, 1 - x - y - z};
		  for (int j = 0; j < 4; j++)
		  {
			  Point < 3> p = mesh[el[j]];
			  for (int k = 0; k < 3; k++)
			  {
				mapped_x[k] += p(k) * lami[j];
			  }
		  }
		  if (info.order == 1)
		  {
			  break;
		  }

	  ELEMENT_EDGE edges = MeshTopology.GetEdges1(TET);
	  for (int i = 0; i < 6; i++)
	  {
		  int eorder = edgeorder[info.edgenrs[i]];
		  if (eorder >= 2)
		  {
				  int first = edgecoeffsindex[info.edgenrs[i]];

		  int vi1 = edges[i][0] - 1;
		  int vi2 = edges[i][1] - 1;
		  if (el[vi1] > el[vi2])
		  {
			  swap(vi1, vi2);
		  }

		  netgen.GlobalMembers.CalcScaledEdgeShapeLambda(eorder, lami[vi1] - lami[vi2], lami[vi1] + lami[vi2], (int i, AutoDiff<3,T> shape) =>
		  {
											   Vec < 3> coef = edgecoeffs[first + i];
											   for (int k = 0; k < 3; k++)
											   {
												 mapped_x[k] += coef(k) * shape;
											   }
		  });
		  }
	  }

	  ELEMENT_FACE faces = MeshTopology.GetFaces1(TET);
	  for (int i = 0; i < 4; i++)
	  {
		  int forder = faceorder[info.facenrs[i]];
		  if (forder >= 3)
		  {
				  int first = facecoeffsindex[info.facenrs[i]];

		  int[] fnums = {faces[i][0] - 1, faces[i][1] - 1, faces[i][2] - 1};
		  if (el[fnums[0]] > el[fnums[1]])
		  {
			  swap(fnums[0], fnums[1]);
		  }
		  if (el[fnums[1]] > el[fnums[2]])
		  {
			  swap(fnums[1], fnums[2]);
		  }
		  if (el[fnums[0]] > el[fnums[1]])
		  {
			  swap(fnums[0], fnums[1]);
		  }

		  netgen.GlobalMembers.CalcScaledTrigShapeLambda(forder, lami[fnums[1]] - lami[fnums[0]], lami[fnums[2]], lami[fnums[0]] + lami[fnums[1]] + lami[fnums[2]], (int i, AutoDiff<3,T> shape) =>
		  {
											   Vec < 3> coef = facecoeffs[first + i];
											   for (int k = 0; k < 3; k++)
											   {
												 mapped_x[k] += coef(k) * shape;
											   }
		  });
		  }
	  }

		  break;
	  }
	  case HEX:
	  {
		  if (info.order >= 2)
		  {
			  return false; // not yet supported
		  }
//C++ TO C# CONVERTER TODO TASK: The following line could not be converted:
		  AutoDiff<3,T> lami[8] = {(1 - x) * (1 - y) * (1 - z), (x) * (1 - y) * (1 - z), (x) * y * (1 - z), (1 - x) * y * (1 - z), (1 - x) * (1 - y) * (z), (x) * (1 - y) * (z), (x) * y * (z), (1 - x) * y * (z)};

		  for (int j = 0; j < 8; j++)
		  {
			  Point < 3> p = mesh[el[j]];
			  for (int k = 0; k < 3; k++)
			  {
				mapped_x[k] += p(k) * lami[j];
			  }
		  }

	  if (info.order == 1)
	  {
		  break;
	  }

//C++ TO C# CONVERTER TODO TASK: The following line could not be converted:
	  AutoDiff<3,T> mu[8] = {(1 - x) + (1 - y) + (1 - z), x + (1 - y) + (1 - z), x + y + (1 - z), (1 - x) + y + (1 - z), (1 - x) + (1 - y) + (z), x + (1 - y) + (z), x + y + (z), (1 - x) + y + (z)};
	  // int ii = 8;
	  ELEMENT_EDGE edges = MeshTopology.GetEdges1(HEX);

	  for (int i = 0; i < 8; i++)
	  {
		  int eorder = edgeorder[info.edgenrs[i]];
		  if (eorder >= 2)
		  {
				  int first = edgecoeffsindex[info.edgenrs[i]];
		  int vi1 = edges[i][0] - 1;
		  int vi2 = edges[i][1] - 1;
		  if (el[vi1] > el[vi2])
		  {
			  swap(vi1, vi2);
		  }

				  AutoDiff<3,T> lame = lami[vi1] + lami[vi2];
		  netgen.GlobalMembers.CalcEdgeShapeLambda(eorder, mu[vi1] - mu[vi2], (int i, AutoDiff<3,T> shape) =>
		  {
										 Vec < 3> coef = edgecoeffs[first + i];
										 for (int k = 0; k < 3; k++)
										 {
										   mapped_x[k] += coef(k) * (lame * shape);
										 }
		  });

		  }
	  }

		  break;
	  }
	  default:
		return false;
	}

	for (int i = 0; i < 3; i++)
	{
		mx(i) = mapped_x[i].Value();
		for (int j = 0; j < 3; j++)
		{
		  jac(i,j) = mapped_x[i].DValue(j);
		}
	}
	return true;
  }

  private class SurfaceElementInfo
  {
	public SurfaceElementIndex elnr = new SurfaceElementIndex();
	public int order;
	public int nv;
	public int ndof;
	public ArrayMem<int,4> edgenrs = new ArrayMem<int,4>();
	public int facenr;
  }

//C++ TO C# CONVERTER TODO TASK: The original C++ template specifier was replaced with a C# generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <typename T>
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void CalcElementShapes(SurfaceElementInfo & info, const Point<2,T> xi, TFlatVector<T> shapes) const
  private void CalcElementShapes<T>(SurfaceElementInfo info, Point<2,T> xi, TFlatVector<T> shapes)
  {
	Element2d el = mesh[info.elnr];
	// shapes.SetSize(info.ndof);

	if (rational && info.order >= 2)
	{
	// shapes.SetSize(6);
	T w = new T(1);
	T[] lami = {xi(0), xi(1), 1 - xi(0) - xi(1)};
	for (int j = 0; j < 3; j++)
	{
	  shapes(j) = lami[j] * lami[j];
	}

	ELEMENT_EDGE edges = MeshTopology.GetEdges1(TRIG);
	for (int j = 0; j < 3; j++)
	{
		T wi = edgeweight[info.edgenrs[j]];
		shapes(j + 3) = 2 * wi * lami[edges[j][0] - 1] * lami[edges[j][1] - 1];
		w += (wi - 1) * 2 * lami[edges[j][0] - 1] * lami[edges[j][1] - 1];
	}

	shapes *= 1.0 / w;
	return;
	}

	switch (el.GetType())
	{
	  case TRIG:
	  {
	  shapes(0) = xi(0);
	  shapes(1) = xi(1);
	  shapes(2) = 1 - xi(0) - xi(1);

	  if (info.order == 1)
	  {
		  return;
	  }

	  int ii = 3;
	  ELEMENT_EDGE edges = MeshTopology.GetEdges0(TRIG);

	  for (int i = 0; i < 3; i++)
	  {
		  int eorder = edgeorder[info.edgenrs[i]];
		  if (eorder >= 2)
		  {
		  int vi1 = edges[i][0];
		  int vi2 = edges[i][1];
		  if (el[vi1] > el[vi2])
		  {
			  swap(vi1, vi2);
		  }

		  netgen.GlobalMembers.CalcScaledEdgeShape(eorder, shapes(vi1) - shapes(vi2), shapes(vi1) + shapes(vi2), shapes(ii));
		  ii += eorder - 1;
		  }
	  }

	  int forder = faceorder[info.facenr];
	  if (forder >= 3)
	  {
		  int[] fnums = {0, 1, 2};
		  if (el[fnums[0]] > el[fnums[1]])
		  {
			  swap(fnums[0], fnums[1]);
		  }
		  if (el[fnums[1]] > el[fnums[2]])
		  {
			  swap(fnums[1], fnums[2]);
		  }
		  if (el[fnums[0]] > el[fnums[1]])
		  {
			  swap(fnums[0], fnums[1]);
		  }

		  netgen.GlobalMembers.CalcTrigShape(forder, shapes(fnums[1]) - shapes(fnums[0]), 1 - shapes(fnums[1]) - shapes(fnums[0]), shapes(ii));
	  }
	  break;
	  }

	  case TRIG6:
	  {
	  if (shapes.Size() == 3)
	  {
		  shapes(0) = xi(0);
		  shapes(1) = xi(1);
		  shapes(2) = 1 - xi(0) - xi(1);
	  }
	  else
	  {
		  T x = xi(0);
		  T y = xi(1);
		  T lam3 = 1 - x - y;

		  shapes(0) = x * (2 * x - 1);
		  shapes(1) = y * (2 * y - 1);
		  shapes(2) = lam3 * (2 * lam3 - 1);
		  shapes(3) = 4 * y * lam3;
		  shapes(4) = 4 * x * lam3;
		  shapes(5) = 4 * x * y;
	  }
	  break;
	  }

	  case QUAD:
	  {
	  shapes(0) = (1 - xi(0)) * (1 - xi(1));
	  shapes(1) = xi(0) * (1 - xi(1));
	  shapes(2) = xi(0) * xi(1);
	  shapes(3) = (1 - xi(0)) * xi(1);

	  if (info.order == 1)
	  {
		  return;
	  }

	  T[] mu = {1 - xi(0) + 1 - xi(1), xi(0) + 1 - xi(1), xi(0) + xi(1), 1 - xi(0) + xi(1)};

	  int ii = 4;
	  ELEMENT_EDGE edges = MeshTopology.GetEdges1(QUAD);

	  for (int i = 0; i < 4; i++)
	  {
		  int eorder = edgeorder[info.edgenrs[i]];
		  if (eorder >= 2)
		  {
		  int vi1 = edges[i][0] - 1;
		  int vi2 = edges[i][1] - 1;
		  if (el[vi1] > el[vi2])
		  {
			  swap(vi1, vi2);
		  }

		  netgen.GlobalMembers.CalcEdgeShape(eorder, mu[vi1] - mu[vi2], shapes(ii));
		  T lame = shapes(vi1) + shapes(vi2);
		  for (int j = 0; j < order - 1; j++)
		  {
			shapes(ii + j) *= lame;
		  }
		  ii += eorder - 1;
		  }
	  }

	  for (int i = ii; i < info.ndof; i++)
	  {
		shapes(i) = 0;
	  }

	  break;
	  }

	  case QUAD8:
	  {
		  var x = xi(0);
		  var y = xi(1);
	  shapes(0) = (1 - x) * (1 - y);
	  shapes(1) = x * (1 - y);
	  shapes(2) = x * y;
	  shapes(3) = (1 - x) * y;
		  shapes(4) = 4 * (1 - x) * x * (1 - y);
		  shapes(5) = 4 * (1 - x) * x * y;
		  shapes(6) = 4 * (1 - y) * y * (1 - x);
		  shapes(7) = 4 * (1 - y) * y * x;
		  shapes(0) -= 0.5 * (shapes(4) + shapes(6));
		  shapes(1) -= 0.5 * (shapes(4) + shapes(7));
		  shapes(2) -= 0.5 * (shapes(5) + shapes(7));
		  shapes(3) -= 0.5 * (shapes(5) + shapes(6));
		  break;
	  }

	  default:
	throw new Exception("CurvedElements::CalcShape 2d, element type not handled");
	};
  }

//C++ TO C# CONVERTER TODO TASK: C++ template specifiers with non-type parameters cannot be converted to C#:
//ORIGINAL LINE: template <int DIM_SPACE>
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetCoefficients(SurfaceElementInfo & info, Array<Vec<DIM_SPACE>> & coefs) const
  private void GetCoefficients<int DIM_SPACE>(SurfaceElementInfo info, Array<Vec<DIM_SPACE>> coefs)
  {
	Element2d el = mesh[info.elnr];
	coefs.SetSize(info.ndof);

	for (int i = 0; i < info.nv; i++)
	{
	Point < 3> hv = mesh[el[i]];
	for (int j = 0; j < DIM_SPACE; j++)
	{
	  coefs[i](j) = hv(j);
	}
	}

	if (info.order == 1)
	{
		return;
	}

	int ii = info.nv;

	for (int i = 0; i < info.edgenrs.Size(); i++)
	{
	int first = edgecoeffsindex[info.edgenrs[i]];
	int next = edgecoeffsindex[info.edgenrs[i] + 1];
	for (int j = first; j < next; j++, ii++)
	{
	  for (int k = 0; k < DIM_SPACE; k++)
	  {
		coefs[ii](k) = edgecoeffs[j](k);
	  }
	}
	}

	int first = facecoeffsindex[info.facenr];
	int next = facecoeffsindex[info.facenr + 1];
	for (int j = first; j < next; j++, ii++)
	{
	  for (int k = 0; k < DIM_SPACE; k++)
	  {
	coefs[ii](k) = facecoeffs[j](k);
	  }
	}
  }

//C++ TO C# CONVERTER TODO TASK: The original C++ template specifier was replaced with a C# generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <typename T>
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void CalcElementDShapes(SurfaceElementInfo & info, const Point<2,T> xi, MatrixFixWidth<2,T> dshapes) const
  private void CalcElementDShapes<T>(SurfaceElementInfo info, Point<2,T> xi, MatrixFixWidth<2,T> dshapes)
  {
	Element2d el = mesh[info.elnr];
	ELEMENT_TYPE type = el.GetType();

	T[] lami = Arrays.InitializeWithDefaultInstances<T>(4);

	dshapes.SetSize(info.ndof);
	// dshapes = 0;

	// *testout << "calcelementdshapes, info.ndof = " << info.ndof << endl;

	if (rational && info.order >= 2)
	{
	T w = 1;
	T[] dw = {0, 0};


	lami[0] = xi(0);
	lami[1] = xi(1);
	lami[2] = 1 - xi(0) - xi(1);
	T[][] dlami =
	{
		new T[] {1, 0},
		new T[] {0, 1},
		new T[] {-1, -1}
	};
	T[] shapes = Arrays.InitializeWithDefaultInstances<T>(6);

	for (int j = 0; j < 3; j++)
	{
		shapes[j] = lami[j] * lami[j];
		dshapes(j,0) = 2 * lami[j] * dlami[j][0];
		dshapes(j,1) = 2 * lami[j] * dlami[j][1];
	}

	ELEMENT_EDGE edges = MeshTopology.GetEdges1(TRIG);
	for (int j = 0; j < 3; j++)
	{
		T wi = edgeweight[info.edgenrs[j]];

		shapes[j + 3] = 2 * wi * lami[edges[j][0] - 1] * lami[edges[j][1] - 1];
		for (int k = 0; k < 2; k++)
		{
		  dshapes(j + 3,k) = 2 * wi * (lami[edges[j][0] - 1] * dlami[edges[j][1] - 1][k] + lami[edges[j][1] - 1] * dlami[edges[j][0] - 1][k]);
		}

		w += (wi - 1) * 2 * lami[edges[j][0] - 1] * lami[edges[j][1] - 1];
		for (int k = 0; k < 2; k++)
		{
		  dw[k] += 2 * (wi - 1) * (lami[edges[j][0] - 1] * dlami[edges[j][1] - 1][k] + lami[edges[j][1] - 1] * dlami[edges[j][0] - 1][k]);
		}
	}
	// shapes *= 1.0 / w;
	dshapes *= 1.0 / w;
	for (int i = 0; i < 6; i++)
	{
	  for (int j = 0; j < 2; j++)
	  {
		dshapes(i,j) -= shapes[i] * dw[j] / (w * w);
	  }
	}
	return;
	}





	switch (type)
	{
	  case TRIG:
	  {
	  dshapes(0,0) = 1;
	  dshapes(0,1) = 0.0;
	  dshapes(1,0) = 0.0;
	  dshapes(1,1) = 1;
	  dshapes(2,0) = -1;
	  dshapes(2,1) = -1;

	  if (info.order == 1)
	  {
		  return;
	  }

	  // *testout << "info.order = " << info.order << endl;


	  lami[0] = xi(0);
	  lami[1] = xi(1);
	  lami[2] = 1 - xi(0) - xi(1);

	  int ii = 3;
	  ELEMENT_EDGE edges = MeshTopology.GetEdges1(TRIG);

	  for (int i = 0; i < 3; i++)
	  {
		  int eorder = edgeorder[info.edgenrs[i]];
		  if (eorder >= 2)
		  {
		  int vi1 = edges[i][0] - 1;
		  int vi2 = edges[i][1] - 1;
		  if (el[vi1] > el[vi2])
		  {
			  swap(vi1, vi2);
		  }

		  netgen.GlobalMembers.CalcScaledEdgeShapeDxDt < 2> (eorder, lami[vi1] - lami[vi2], lami[vi1] + lami[vi2], dshapes(ii,0));

		  Mat<2,2,T> trans = new Mat<2,2,T>();
		  for (int j = 0; j < 2; j++)
		  {
			  trans(0,j) = dshapes(vi1,j) - dshapes(vi2,j);
			  trans(1,j) = dshapes(vi1,j) + dshapes(vi2,j);
		  }

		  for (int j = 0; j < eorder - 1; j++)
		  {
			  T ddx = dshapes(ii + j,0);
			  T ddt = dshapes(ii + j,1);
			  dshapes(ii + j,0) = ddx * trans(0,0) + ddt * trans(1,0);
			  dshapes(ii + j,1) = ddx * trans(0,1) + ddt * trans(1,1);
		  }

		  ii += eorder - 1;
		  }
	  }

	  int forder = faceorder[info.facenr];
	  // *testout << "forder = " << forder << endl;
	  if (forder >= 3)
	  {
		  int[] fnums = {0, 1, 2};
		  if (el[fnums[0]] > el[fnums[1]])
		  {
			  swap(fnums[0], fnums[1]);
		  }
		  if (el[fnums[1]] > el[fnums[2]])
		  {
			  swap(fnums[1], fnums[2]);
		  }
		  if (el[fnums[0]] > el[fnums[1]])
		  {
			  swap(fnums[0], fnums[1]);
		  }

		  netgen.GlobalMembers.CalcTrigShapeDxDy(forder, lami[fnums[1]] - lami[fnums[0]], 1 - lami[fnums[1]] - lami[fnums[0]], dshapes(ii,0));

		  int nd = (forder - 1) * (forder - 2) / 2;
		  Mat<2,2,T> trans = new Mat<2,2,T>();
		  for (int j = 0; j < 2; j++)
		  {
		  trans(0,j) = dshapes(fnums[1],j) - dshapes(fnums[0],j);
		  trans(1,j) = -dshapes(fnums[1],j) - dshapes(fnums[0],j);
		  }

		  for (int j = 0; j < nd; j++)
		  {
		  T ddx = dshapes(ii + j,0);
		  T ddt = dshapes(ii + j,1);
		  dshapes(ii + j,0) = ddx * trans(0,0) + ddt * trans(1,0);
		  dshapes(ii + j,1) = ddx * trans(0,1) + ddt * trans(1,1);
		  }
	  }

	  break;
	  }

	  case TRIG6:
	  {
	  if (dshapes.Height() == 3)
	  {
		  dshapes = T(0.0);
		  dshapes(0,0) = 1;
		  dshapes(1,1) = 1;
		  dshapes(2,0) = -1;
		  dshapes(2,1) = -1;
	  }
	  else
	  {
		  AutoDiff<2,T> x = new AutoDiff<2,T>(xi(0), 0);
		  AutoDiff<2,T> y = new AutoDiff<2,T>(xi(1), 1);
		  AutoDiff<2,T> lam3 = 1 - x - y;
		  AutoDiff<2,T>[] shapes = Arrays.InitializeWithDefaultInstances<AutoDiff>(6);
		  shapes[0] = x * (2 * x - 1);
		  shapes[1] = y * (2 * y - 1);
		  shapes[2] = lam3 * (2 * lam3 - 1);
		  shapes[3] = 4 * y * lam3;
		  shapes[4] = 4 * x * lam3;
		  shapes[5] = 4 * x * y;

		  for (int i = 0; i < 6; i++)
		  {
		  dshapes(i,0) = shapes[i].DValue(0);
		  dshapes(i,1) = shapes[i].DValue(1);
		  }

	  }
	  break;
	  }

	  case QUAD:
	  {
	  dshapes(0,0) = -(1 - xi(1));
	  dshapes(0,1) = -(1 - xi(0));
	  dshapes(1,0) = (1 - xi(1));
	  dshapes(1,1) = -xi(0);
	  dshapes(2,0) = xi(1);
	  dshapes(2,1) = xi(0);
	  dshapes(3,0) = -xi(1);
	  dshapes(3,1) = (1 - xi(0));

	  if (info.order == 1)
	  {
		  return;
	  }

	  T[] shapes = {(1 - xi(0)) * (1 - xi(1)), xi(0) * (1 - xi(1)), xi(0) * xi(1), (1 - xi(0)) * xi(1)};

	  T[] mu = {1 - xi(0) + 1 - xi(1), xi(0) + 1 - xi(1), xi(0) + xi(1), 1 - xi(0) + xi(1)};

	  T[][] dmu =
	  {
		  new T[] {-1, -1},
		  new T[] {1, -1},
		  new T[] {1, 1},
		  new T[] {-1, 1}
	  };

	  // double hshapes[20], hdshapes[20];
	  ArrayMem<T, 20> hshapes = new ArrayMem<T, 20>(order + 1);
	  ArrayMem<T, 20> hdshapes = new ArrayMem<T, 20>(order + 1);

	  int ii = 4;
	  ELEMENT_EDGE edges = MeshTopology.GetEdges1(QUAD);

	  for (int i = 0; i < 4; i++)
	  {
		  int eorder = edgeorder[info.edgenrs[i]];
		  if (eorder >= 2)
		  {
		  int vi1 = edges[i][0] - 1;
		  int vi2 = edges[i][1] - 1;
		  if (el[vi1] > el[vi2])
		  {
			  swap(vi1, vi2);
		  }

		  netgen.GlobalMembers.CalcEdgeShapeDx(eorder, mu[vi1] - mu[vi2], hshapes[0], hdshapes[0]);

		  T lame = shapes[vi1] + shapes[vi2];
		  T[] dlame = {dshapes(vi1, 0) + dshapes(vi2, 0), dshapes(vi1, 1) + dshapes(vi2, 1)};

		  for (int j = 0; j < eorder - 1; j++)
		  {
			for (int k = 0; k < 2; k++)
			{
			  dshapes(ii + j, k) = lame * hdshapes[j] * (dmu[vi1][k] - dmu[vi2][k]) + dlame[k] * hshapes[j];
			}
		  }

		  ii += eorder - 1;
		  }
	  }

	  /*
	   *testout << "quad, dshape = " << endl << dshapes << endl;
	   for (int i = 0; i < 2; i++)
	   {
	   Point<2> xil = xi, xir = xi;
	   Vector shapesl(dshapes.Height()), shapesr(dshapes.Height());
	   xil(i) -= 1e-6;
	   xir(i) += 1e-6;
	   CalcElementShapes (info, xil, shapesl);
	   CalcElementShapes (info, xir, shapesr);
	  
	   for (int j = 0; j < dshapes.Height(); j++)
	   dshapes(j,i) = 1.0 / 2e-6 * (shapesr(j)-shapesl(j));
	   }
	  
	   *testout << "quad, num dshape = " << endl << dshapes << endl;
	   */
	  break;
	  }
	  default:
	throw new Exception("CurvedElements::CalcDShape 2d, element type not handled");

	};
  }

//C++ TO C# CONVERTER TODO TASK: C++ template specifiers with non-type parameters cannot be converted to C#:
//ORIGINAL LINE: template <int DIM_SPACE, typename T>
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool EvaluateMapping(SurfaceElementInfo & info, const Point<2,T> xi, Point<DIM_SPACE,T> & mx, Mat<DIM_SPACE,2,T> & jac) const
  private bool EvaluateMapping<int DIM_SPACE, T>(SurfaceElementInfo info, Point<2,T> xi, Point<DIM_SPACE,T> mx, Mat<DIM_SPACE,2,T> jac)
  {
	Element2d el = mesh[info.elnr];
	if (rational && info.order >= 2)
	{
		return false; // not supported
	}

	AutoDiff<2,T> x = new AutoDiff<2,T>(xi(0), 0);
	AutoDiff<2,T> y = new AutoDiff<2,T>(xi(1), 1);

	AutoDiff<2,T>[] mapped_x = Arrays.InitializeWithDefaultInstances<AutoDiff>(DIM_SPACE);
	for (int i = 0; i < DIM_SPACE; i++)
	{
	  mapped_x[i] = new AutoDiff<2,T>(0.0);
	}

	switch (el.GetType())
	{
	  case TRIG:
	  {
		  // if (info.order >= 2) return false; // not yet supported
//C++ TO C# CONVERTER TODO TASK: The following line could not be converted:
		  AutoDiff<2,T> lami[4] = {x, y, 1 - x - y};
		  for (int j = 0; j < 3; j++)
		  {
			  Point < 3> p = mesh[el[j]];
			  for (int k = 0; k < DIM_SPACE; k++)
			  {
				mapped_x[k] += p(k) * lami[j];
			  }
		  }
		  if (info.order == 1)
		  {
			  break;
		  }

	  ELEMENT_EDGE edges = MeshTopology.GetEdges1(TRIG);
	  for (int i = 0; i < 3; i++)
	  {
		  int eorder = edgeorder[info.edgenrs[i]];
		  if (eorder >= 2)
		  {
				  int first = edgecoeffsindex[info.edgenrs[i]];

		  int vi1 = edges[i][0] - 1;
		  int vi2 = edges[i][1] - 1;
		  if (el[vi1] > el[vi2])
		  {
			  swap(vi1, vi2);
		  }

		  netgen.GlobalMembers.CalcScaledEdgeShapeLambda(eorder, lami[vi1] - lami[vi2], lami[vi1] + lami[vi2], (int i, AutoDiff<2,T> shape) =>
		  {
											   for (int k = 0; k < DIM_SPACE; k++)
											   {
												 mapped_x[k] += edgecoeffs[first + i](k) * shape;
											   }
		  });
		  }
	  }

		  int forder = faceorder[info.facenr];
		  if (forder >= 3)
		  {
			  int first = facecoeffsindex[info.facenr];

			  int[] fnums = {0, 1, 2};
			  if (el[fnums[0]] > el[fnums[1]])
			  {
				  swap(fnums[0], fnums[1]);
			  }
			  if (el[fnums[1]] > el[fnums[2]])
			  {
				  swap(fnums[1], fnums[2]);
			  }
			  if (el[fnums[0]] > el[fnums[1]])
			  {
				  swap(fnums[0], fnums[1]);
			  }

			  netgen.GlobalMembers.CalcScaledTrigShapeLambda(forder, lami[fnums[1]] - lami[fnums[0]], lami[fnums[2]], new AutoDiff<2,T>(1.0), (int i, AutoDiff<2,T> shape) =>
			  {
										   for (int k = 0; k < DIM_SPACE; k++)
										   {
											 mapped_x[k] += facecoeffs[first + i](k) * shape;
										   }
			  });
		  }
		  break;
	  }
	  case QUAD:
	  {
		  if (info.order >= 2)
		  {
			  return false; // not yet supported
		  }
//C++ TO C# CONVERTER TODO TASK: The following line could not be converted:
		  AutoDiff<2,T> lami[4] = {(1 - x) * (1 - y), x * (1 - y), x * y, (1 - x) * y};
		  for (int j = 0; j < 4; j++)
		  {
			  Point < 3> p = mesh[el[j]];
			  for (int k = 0; k < DIM_SPACE; k++)
			  {
				mapped_x[k] += p(k) * lami[j];
			  }
		  }
		  break;
	  }
	  case QUAD8:
	  {
		  // AutoDiff<2,T> lami[4] = { (1-x)*(1-y), x*(1-y), x*y, (1-x)*y };
//C++ TO C# CONVERTER TODO TASK: The following line could not be converted:
		  AutoDiff<2,T> lami[8] = {(1 - x) * (1 - y), x * (1 - y), x * y, (1 - x) * y, 4 * (1 - x) * x * (1 - y), 4 * (1 - x) * x * y, 4 * (1 - y) * y * (1 - x), 4 * (1 - y) * y * x};

		  lami[0] -= 0.5 * (lami[4] + lami[6]);
		  lami[1] -= 0.5 * (lami[4] + lami[7]);
		  lami[2] -= 0.5 * (lami[5] + lami[7]);
		  lami[3] -= 0.5 * (lami[5] + lami[6]);

		  for (int j = 0; j < 8; j++)
		  {
			  Point < 3> p = mesh[el[j]];
			  for (int k = 0; k < DIM_SPACE; k++)
			  {
				mapped_x[k] += p(k) * lami[j];
			  }
		  }
		  break;
	  }

	  default:
		return false;
	}

	for (int i = 0; i < DIM_SPACE; i++)
	{
		mx(i) = mapped_x[i].Value();
		for (int j = 0; j < 2; j++)
		{
		  jac(i,j) = mapped_x[i].DValue(j);
		}
	}
	return true;
  }
}







namespace netgen
{

  public class RecPol : System.IDisposable
  {
	protected int maxorder;
	protected double[] a;
	protected double[] b;
	protected double[] c;
	public RecPol(int amaxorder)
	{
	  maxorder = amaxorder;
	  a = new double[maxorder + 1];
	  b = new double[maxorder + 1];
	  c = new double[maxorder + 1];
	}
	public void Dispose()
	{
	  Arrays.DeleteArray(a);
	  Arrays.DeleteArray(b);
	  Arrays.DeleteArray(c);
	}

//C++ TO C# CONVERTER TODO TASK: The original C++ template specifier was replaced with a C# generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class S, class T>
	public void Evaluate<S, T>(int n, S x, T[] values)
	{
	  S p1 = new S(1.0);
	  S p2 = new S(0.0);
	  S p3 = new default(S);

	  if (n >= 0)
	  {
	p2 = values[0] = 1.0;
	  }
	  if (n >= 1)
	  {
	p1 = values[1] = a[0] + b[0] * x;
	  }

	  for (int i = 1; i < n; i++)
	  {
	  p3 = p2;
	  p2 = p1;
	  p1 = (a[i] + b[i] * x) * p2 - c[i] * p3;
	  values[i + 1] = p1;
	  }
	}

//C++ TO C# CONVERTER TODO TASK: The original C++ template specifier was replaced with a C# generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class S, class T>
	public void EvaluateScaled<S, T>(int n, S x, S y, T[] values)
	{
	  S p1 = new S(1.0);
	  S p2 = new S(0.0);
	  S p3 = new default(S);

	  if (n >= 0)
	  {
	p2 = values[0] = 1.0;
	  }
	  if (n >= 1)
	  {
	p1 = values[1] = a[0] * y + b[0] * x;
	  }

	  for (int i = 1; i < n; i++)
	  {
	  p3 = p2;
	  p2 = p1;
	  p1 = (a[i] * y + b[i] * x) * p2 - c[i] * y * y * p3;
	  values[i + 1] = p1;
	  }
	}

//C++ TO C# CONVERTER TODO TASK: The original C++ template specifier was replaced with a C# generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class S, class FUNC>
	public void EvaluateScaledLambda<S, FUNC>(int n, S x, S y, FUNC func)
	{
	  S p1 = new S(1.0);
	  S p2 = new S(0.0);
	  S p3 = new default(S);

	  if (n >= 0)
	  {
		  p2 = 1.0;
		  func(0, p2);
	  }
	  if (n >= 1)
	  {
		  p1 = a[0] * y + b[0] * x;
		  func(1, p1);
	  }

	  for (int i = 1; i < n; i++)
	  {
	  p3 = p2;
	  p2 = p1;
	  p1 = (a[i] * y + b[i] * x) * p2 - c[i] * y * y * p3;
	  func(i + 1, p1);
	  }
	}


  }

  public class JacobiRecPol : RecPol
  {
	public JacobiRecPol(int amo, double al, double be) : base(amo)
	{
	  for (int i = 0; i <= maxorder; i++)
	  {
	  double den = 2 * (i + 1) * (i + al + be+1) * (2 * i + al + be);
	  a[i] = (2 * i + al + be+1) * (al * al - be * be) / den;
	  b[i] = (2 * i + al + be) * (2 * i + al + be+1) * (2 * i + al + be+2) / den;
	  c[i] = 2 * (i + al) * (i + be) * (2 * i + al + be+2) / den;
	  }
	}
  }



//C++ TO C# CONVERTER TODO TASK: The original C++ template specifier was replaced with a C# generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class S, class T>

}
