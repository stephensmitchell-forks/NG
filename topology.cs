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
/* File:   topology.hh                                                    */
/* Author: Joachim Schoeberl                                              */
/* Date:   27. Apr. 01                                                    */
/**************************************************************************/

/*
    Mesh topology
    (Elements, Faces, Edges, Vertices
*/

namespace netgen
{

public class T_EDGE
{
  // int orient:1;
  public int nr; // 0-based
}

public class T_FACE
{
  // int forient:3;
  public int fnr; // 0-based
}

//C++ TO C# CONVERTER TODO TASK: C++ template specifiers with non-type parameters cannot be converted to C#:
//ORIGINAL LINE: template <typename T, int S>
  public class FixArray <T, int S>
  {
	public T[] vals = Arrays.InitializeWithDefaultInstances<T>(S);
	public T this [uint i]
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
//ORIGINAL LINE: T operator [] (uint i) const
	public T this [uint i]
	{
		get
		{
			return vals[i];
		}
	}
  }


public class MeshTopology : System.IDisposable
{
  private readonly Mesh[] mesh;
  private bool buildedges;
  private bool buildfaces;

  private Array<INDEX_2> edge2vert = new Array<INDEX_2>();
  private Array<INDEX_4> face2vert = new Array<INDEX_4>();
  /*
  Array<T_EDGE[12]> edges;
  Array<T_FACE[6]> faces;
  Array<T_EDGE[4]> surfedges;
  */
  private Array<FixArray<T_EDGE,12>> edges = new Array<FixArray<T_EDGE,12>>();
  private Array<FixArray<T_FACE,6>> faces = new Array<FixArray<T_FACE,6>>();
  private Array<FixArray<T_EDGE,4>> surfedges = new Array<FixArray<T_EDGE,4>>();

  private Array<T_EDGE> segedges = new Array<T_EDGE>();
  private Array<T_FACE> surffaces = new Array<T_FACE>();
  private Array<INDEX_2> surf2volelement = new Array<INDEX_2>();
  private Array<int> face2surfel = new Array<int>();
  private TABLE<ElementIndex,PointIndex.BASE> vert2element = new TABLE<ElementIndex,PointIndex.BASE>();
  private TABLE<SurfaceElementIndex,PointIndex.BASE> vert2surfelement = new TABLE<SurfaceElementIndex,PointIndex.BASE>();
  private TABLE<SegmentIndex,PointIndex.BASE> vert2segment = new TABLE<SegmentIndex,PointIndex.BASE>();
  private TABLE<int,PointIndex.BASE> vert2pointelement = new TABLE<int,PointIndex.BASE>();
  private int timestamp;
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetNSurfedges() const
  public int GetNSurfedges()
  {
	  return surfedges.Size();
  }

//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//  MeshTopology() = default;
//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//  MeshTopology(const MeshTopology & top) = default;
//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//  MeshTopology(MeshTopology && top) = default;
  public MeshTopology(Mesh amesh)
  {
	  this.mesh = new netgen.Mesh(amesh);
	buildedges = true;
	buildfaces = true;
	timestamp = -1;
  }

  public void Dispose()
  {
	  ;
  }

//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//  MeshTopology & operator = (const MeshTopology & top) = default;
//C++ TO C# CONVERTER TODO TASK: C# has no equivalent to ' = default':
//  MeshTopology & operator = (MeshTopology && top) = default;

  public void SetBuildEdges(bool be)
  {
	  buildedges = be;
  }
  public void SetBuildFaces(bool bf)
  {
	  buildfaces = bf;
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool HasEdges() const
  public bool HasEdges()
  {
	  return buildedges;
  }
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool HasFaces() const
  public bool HasFaces()
  {
	  return buildfaces;
  }

//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private int Update_timer = NgProfiler.CreateTimer("topology");
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private int Update_timer1 = NgProfiler.CreateTimer("topology::buildedges");
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private int Update_timer2 = NgProfiler.CreateTimer("topology::buildfaces");
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private int Update_timer2a = NgProfiler.CreateTimer("topology::buildfacesa");
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private int Update_timer2b = NgProfiler.CreateTimer("topology::buildfacesb");
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private int Update_timer2b1 = NgProfiler.CreateTimer("topology::buildfacesb1");
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private int Update_timer2c = NgProfiler.CreateTimer("topology::buildfacesc");

  public void Update(TaskManager tm = DummyTaskManager, Tracer tracer = DummyTracer)
  {
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static int timer = NgProfiler::CreateTimer("topology");
	NgProfiler.RegionTimer reg = new NgProfiler.RegionTimer(Update_timer);

#if PARALLEL
	// ParallelMeshTopology & paralleltop = mesh.GetParallelTopology();
#endif

	var id = this.mesh.GetCommunicator().Rank();
	var ntasks = this.mesh.GetCommunicator().Size();

	if (timestamp > mesh.GetTimeStamp())
	{
		return;
	}

	int ne = mesh.GetNE();
	int nse = mesh.GetNSE();
	int nseg = mesh.GetNSeg();
	int np = mesh.GetNP();
	int nv = mesh.GetNV();

	if (id == 0)
	{
	  PrintMessage(3, "Update mesh topology");
	}

	(*testout) << " UPDATE MESH TOPOLOGY " << "\n";
	(*testout) << "ne   = " << ne << "\n";
	(*testout) << "nse  = " << nse << "\n";
	(*testout) << "nseg = " << nseg << "\n";
	(*testout) << "np   = " << np << "\n";
	(*testout) << "nv   = " << nv << "\n";

	tracer("Topology::Update setup tables", false);
	Array<int,PointIndex.BASE> cnt = new Array<int,PointIndex.BASE>(nv);
	Array<int> vnums = new Array<int>();

	/*
	  generate:
	  vertex to element
	  vertex to surface element
	  vertex to segment
	*/
	cnt = 0;
	/*
	for (ElementIndex ei = 0; ei < ne; ei++)
	  {
	const Element & el = (*mesh)[ei];
	for (int j = 0; j < el.GetNV(); j++)
	  cnt[el[j]]++;
	  }
	*/
	netgen.GlobalMembers.ParallelForRange(tm, (uint)ne, (uint begin, uint end) =>
	{
		 for (ElementIndex ei = begin; ei < end; ei++)
		 {
			 Element el = mesh[ei];
			 for (int j = 0; j < el.GetNV(); j++)
			 {
			   AsAtomic(cnt[el[j]])++;
			 }
		 }
	});

	vert2element = new TABLE<ElementIndex,PointIndex.BASE> (cnt);
	/*
	for (ElementIndex ei = 0; ei < ne; ei++)
	  {
	const Element & el = (*mesh)[ei];
	for (int j = 0; j < el.GetNV(); j++)
	  vert2element.AddSave (el[j], ei);
	  }
	*/
	netgen.GlobalMembers.ParallelForRange(tm, (uint)ne, (uint begin, uint end) =>
	{
		 for (ElementIndex ei = begin; ei < end; ei++)
		 {
			 Element el = mesh[ei];
			 for (int j = 0; j < el.GetNV(); j++)
			 {
			   vert2element.ParallelAdd(el[j], ei);
			 }
		 }
	});

	cnt = 0;
	/*
	for (SurfaceElementIndex sei = 0; sei < nse; sei++)
	  {
	const Element2d & el = (*mesh)[sei];
	for (int j = 0; j < el.GetNV(); j++)
	  cnt[el[j]]++;
	  }
	*/
	netgen.GlobalMembers.ParallelForRange(tm, (uint)nse, (uint begin, uint end) =>
	{
		 for (SurfaceElementIndex ei = begin; ei < end; ei++)
		 {
			 Element2d el = mesh[ei];
			 for (int j = 0; j < el.GetNV(); j++)
			 {
			   AsAtomic(cnt[el[j]])++;
			 }
		 }
	});



	vert2surfelement = new TABLE<SurfaceElementIndex,PointIndex.BASE> (cnt);
	/*
	for (SurfaceElementIndex sei = 0; sei < nse; sei++)
	  {
	const Element2d & el = (*mesh)[sei];
	for (int j = 0; j < el.GetNV(); j++)
	  vert2surfelement.AddSave (el[j], sei);
	  }
	*/
	netgen.GlobalMembers.ParallelForRange(tm, (uint)nse, (uint begin, uint end) =>
	{
		 for (SurfaceElementIndex sei = begin; sei < end; sei++)
		 {
			 Element2d el = mesh[sei];
			 for (int j = 0; j < el.GetNV(); j++)
			 {
			   vert2surfelement.ParallelAdd(el[j], sei);
			 }
		 }
	});


	cnt = 0;
	for (SegmentIndex si = 0; si < nseg; si++)
	{
	Segment seg = mesh.LineSegment(si);
	cnt[seg[0]]++;
	cnt[seg[1]]++;
	}

	vert2segment = new TABLE<SegmentIndex,PointIndex.BASE> (cnt);
	for (SegmentIndex si = 0; si < nseg; si++)
	{
	Segment seg = mesh.LineSegment(si);
	vert2segment.AddSave(seg[0], si);
	vert2segment.AddSave(seg[1], si);
	}


	cnt = 0;
	for (int pei = 0; pei < mesh.pointelements.Size(); pei++)
	{
	Element0d pointel = mesh.pointelements[pei];
	cnt[pointel.pnum]++;
	}

	vert2pointelement = new TABLE<int,PointIndex.BASE> (cnt);
	for (int pei = 0; pei < mesh.pointelements.Size(); pei++)
	{
	Element0d pointel = mesh.pointelements[pei];
	vert2pointelement.AddSave(pointel.pnum, pei);
	}
	tracer("Topology::Update setup tables", true);


	if (buildedges)
	{
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static int timer1 = NgProfiler::CreateTimer("topology::buildedges");
	NgProfiler.RegionTimer reg1 = new NgProfiler.RegionTimer(Update_timer1);

	if (id == 0)
	{
	  PrintMessage(5, "Update edges ");
	}

	edges.SetSize(ne);
	surfedges.SetSize(nse);
	segedges.SetSize(nseg);

	for (int i = 0; i < ne; i++)
	{
	  for (int j = 0; j < 12; j++)
	  {
		edges[i][j].nr = -1;
	  }
	}
	for (int i = 0; i < nse; i++)
	{
	  for (int j = 0; j < 4; j++)
	  {
		surfedges[i][j].nr = -1;
	  }
	}

	// keep existing edges
	cnt = 0;
	for (int i = 0; i < edge2vert.Size(); i++)
	{
	  cnt[edge2vert[i][0]]++;
	}
	TABLE<int,PointIndex.BASE> vert2edge = new TABLE<int,PointIndex.BASE>(cnt);
	for (int i = 0; i < edge2vert.Size(); i++)
	{
	  vert2edge.AddSave(edge2vert[i][0], i);
	}

	// ensure all coarse grid and intermediate level edges
	cnt = 0;
	for (int i = mesh.mlbetweennodes.Begin(); i < mesh.mlbetweennodes.End(); i++)
	{
		INDEX_2 parents = netgen.GlobalMembers.Sort(mesh.mlbetweennodes[i]);
		if (parents[0] >= PointIndex.BASE)
		{
			cnt[parents[0]]++;
		}
	}
	TABLE<int,PointIndex.BASE> vert2vertcoarse = new TABLE<int,PointIndex.BASE>(cnt);
	for (int i = mesh.mlbetweennodes.Begin(); i < mesh.mlbetweennodes.End(); i++)
	{
		INDEX_2 parents = netgen.GlobalMembers.Sort(mesh.mlbetweennodes[i]);
		if (parents[0] >= PointIndex.BASE)
		{
			vert2vertcoarse.AddSave(parents[0], parents[1]);
		}
	}



	int max_edge_on_vertex = 0;
	for (int i = PointIndex.BASE; i < nv + PointIndex.BASE; i++)
	{
		int onv = vert2edge[i].Size() + vert2vertcoarse[i].Size() + 4 * (vert2element)[i].Size() + 2 * (vert2surfelement)[i].Size() + (vert2segment)[i].Size();
		max_edge_on_vertex = Math.Max(onv, max_edge_on_vertex);
	}


		// count edges associated with vertices
		cnt = 0;

		netgen.GlobalMembers.ParallelForRange(tm, (uint)mesh.GetNV(), (uint begin, uint end) =>
		{
			 INDEX_CLOSED_HASHTABLE<int> v2eht = new INDEX_CLOSED_HASHTABLE<int>(2 * max_edge_on_vertex + 10);
			 for (PointIndex v = begin + PointIndex.BASE; v < end + PointIndex.BASE; v++)
			 {
				 v2eht.DeleteData();
				 foreach (int ednr in vert2edge[v])
				 {
					 int v2 = edge2vert[ednr][1];
					 v2eht.Set(v2, ednr);
				 }

				 int cnti = 0;

				 foreach (int v2 in vert2vertcoarse[v])
				 {
				   if (!v2eht.Used(v2))
				   {
					   cnti++;
					   v2eht.Set(v2, 33); // some value
				   }
				 }

				 LoopOverEdges(mesh[0], this, v, (INDEX_2 edge, int elnr, int loc_edge, int element_dim, int edgedir) =>
				 {
								  if (!v2eht.Used(edge.I2()))
								  {
									  cnti++;
									  v2eht.Set(edge.I2(), 33); // something
								  }
				 });
				 cnt[v] = cnti;
			 }
		});

		// accumulate number of edges
		int ned = edge2vert.Size();

		// for (size_t v = 0; v < mesh->GetNV(); v++)
		foreach (uint v in cnt.Range())
		{
			var hv = cnt[v];
			cnt[v] = ned;
			ned += hv;
		}
		edge2vert.SetSize(ned);


		// INDEX_CLOSED_HASHTABLE<int> v2eht(2*max_edge_on_vertex+10);
	// Array<int> vertex2;
	// for (PointIndex v = PointIndex::BASE; v < nv+PointIndex::BASE; v++)

		netgen.GlobalMembers.ParallelForRange(tm, (uint)mesh.GetNV(), (uint begin, uint end) =>
		{
			 INDEX_CLOSED_HASHTABLE<int> v2eht = new INDEX_CLOSED_HASHTABLE<int>(2 * max_edge_on_vertex + 10);
			 Array<int> vertex2 = new Array<int>();
			 for (PointIndex v = begin + PointIndex.BASE; v < end + PointIndex.BASE; v++)
			 {
				 int ned = cnt[v];
				 v2eht.DeleteData();
				 vertex2.SetSize(0);

				 foreach (int ednr in vert2edge[v])
				 {
					 int v2 = edge2vert[ednr][1];
					 v2eht.Set(v2, ednr);
				 }

				 foreach (int v2 in vert2vertcoarse[v])
				 {
				   if (!v2eht.Used(v2))
				   {
					   v2eht.Set(v2, 33); // some value
					   vertex2.Append(v2);
				   }
				 }

				 LoopOverEdges(mesh[0], this, v, (INDEX_2 edge, int elnr, int loc_edge, int element_dim, int edgedir) =>
				 {
								  if (!v2eht.Used(edge.I2()))
								  {
									  vertex2.Append(edge.I2());
									  v2eht.Set(edge.I2(), 33);
								  }
				 });

				 QuickSort(vertex2);

				 for (int j = 0; j < vertex2.Size(); j++)
				 {
					 v2eht.Set(vertex2[j], ned);
					 edge2vert[ned] = new INDEX_2(v, vertex2[j]);
					 ned++;
				 }

				 LoopOverEdges(mesh[0], this, v, (INDEX_2 edge, int elnr, int loc_edge, int element_dim, int edgedir) =>
				 {
								  int edgenum = v2eht.Get(edge.I2());
								  switch (element_dim)
								  {
									case 3:
									  edges[elnr][loc_edge].nr = edgenum;
									  // edges[elnr][loc_edge].orient = edgedir;
									  break;
									case 2:
									  surfedges[elnr][loc_edge].nr = edgenum;
									  // surfedges[elnr][loc_edge].orient = edgedir;
									  break;
									case 1:
									  segedges[elnr].nr = edgenum;
									  // segedges[elnr].orient = edgedir;
									  break;
								  }
				 });
			 }
		});
	}



	// generate faces
	if (buildfaces)
	{
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static int timer2 = NgProfiler::CreateTimer("topology::buildfaces");
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static int timer2a = NgProfiler::CreateTimer("topology::buildfacesa");
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static int timer2b = NgProfiler::CreateTimer("topology::buildfacesb");
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//		static int timer2b1 = NgProfiler::CreateTimer("topology::buildfacesb1");
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static int timer2c = NgProfiler::CreateTimer("topology::buildfacesc");
	NgProfiler.RegionTimer reg2 = new NgProfiler.RegionTimer(Update_timer2);

	if (id == 0)
	{
	  PrintMessage(5, "Update faces ");
	}

		NgProfiler.StartTimer(Update_timer2a);

	faces.SetSize(ne);
	surffaces.SetSize(nse);


	cnt = 0;
	for (int i = 0; i < face2vert.Size(); i++)
	{
	  cnt[face2vert[i][0]]++;
	}
	TABLE<int,PointIndex.BASE> vert2oldface = new TABLE<int,PointIndex.BASE>(cnt);
	for (int i = 0; i < face2vert.Size(); i++)
	{
	  vert2oldface.AddSave(face2vert[i][0], i);
	}


	for (int elnr = 0; elnr < ne; elnr++)
	{
	  for (int j = 0; j < 6; j++)
	  {
		faces[elnr][j].fnr = -1;
	  }
	}


	int max_face_on_vertex = 0;
	for (int i = PointIndex.BASE; i < nv + PointIndex.BASE; i++)
	{
		int onv = vert2oldface[i].Size() + vert2element[i].Size() + vert2surfelement[i].Size();
		max_face_on_vertex = Math.Max(onv, max_face_on_vertex);
	}




		NgProfiler.StopTimer(Update_timer2a);
		NgProfiler.StartTimer(Update_timer2b);

		INDEX_3_CLOSED_HASHTABLE<int> vert2face = new INDEX_3_CLOSED_HASHTABLE<int>(2 * max_face_on_vertex + 10);

	int oldnfa = face2vert.Size();

		// count faces associated with vertices
		cnt = 0;
		// for (auto v : mesh.Points().Range())
		NgProfiler.StartTimer(Update_timer2b1);
		netgen.GlobalMembers.ParallelForRange(tm, (uint)mesh.GetNV(), (uint begin, uint end) =>
		{
			  INDEX_3_CLOSED_HASHTABLE<int> vert2face = new INDEX_3_CLOSED_HASHTABLE<int>(2 * max_face_on_vertex + 10);
			  for (PointIndex v = begin + PointIndex.BASE; v < end + PointIndex.BASE; v++)
			  {
				  vert2face.DeleteData();

				  for (int j = 0; j < vert2oldface[v].Size(); j++)
				  {
					  int fnr = vert2oldface[v][j];
					  INDEX_3 face = new INDEX_3(face2vert[fnr].I1(), face2vert[fnr].I2(), face2vert[fnr].I3());
					  vert2face.Set(face, 33); // something
				  }
				  int cnti = 0;
				  LoopOverFaces(mesh[0], this, v, (INDEX_4 i4, int elnr, int j, bool volume, int facedir) =>
				  {
								   INDEX_3 face = new INDEX_3(i4.I1(), i4.I2(), i4.I3());
								   if (!vert2face.Used(face))
								   {
									   cnti++;
									   vert2face.Set(face, 33); // something
								   }
				  });
				  cnt[v] = cnti;
			  }
		});
		NgProfiler.StopTimer(Update_timer2b1);

		// accumulate number of faces
		int nfa = oldnfa;
		// for (auto v : Range(mesh->GetNV())) // Points().Range())
		// for (size_t v = 0; v < mesh->GetNV(); v++)
		foreach (var v in cnt.Range())
		{
			var hv = cnt[v];
			cnt[v] = nfa;
			nfa += hv;
		}
		face2vert.SetSize(nfa);

		// for (auto v : mesh.Points().Range())

		netgen.GlobalMembers.ParallelForRange(tm, (uint)mesh.GetNV(), (uint begin, uint end) =>
		{
			  INDEX_3_CLOSED_HASHTABLE<int> vert2face = new INDEX_3_CLOSED_HASHTABLE<int>(2 * max_face_on_vertex + 10);
			  for (PointIndex v = begin + PointIndex.BASE; v < end + PointIndex.BASE; v++)
			  {
				  int first_fa = cnt[v];
				  int nfa = first_fa;
				  vert2face.DeleteData();

				  for (int j = 0; j < vert2oldface[v].Size(); j++)
				  {
					  int fnr = vert2oldface[v][j];
					  INDEX_3 face = new INDEX_3(face2vert[fnr].I1(), face2vert[fnr].I2(), face2vert[fnr].I3());
					  vert2face.Set(face, fnr);
				  }

				  LoopOverFaces(mesh[0], this, v, (INDEX_4 i4, int elnr, int j, bool volume, int facedir) =>
				  {
								   INDEX_3 face = new INDEX_3(i4.I1(), i4.I2(), i4.I3());
								   if (!vert2face.Used(face))
								   {
									   face2vert[nfa] = i4;
									   vert2face.Set(face, nfa);
									   nfa++;
								   }
				  });


				  QuickSort(face2vert.Range(first_fa, nfa));

				  for (int j = first_fa; j < nfa; j++)
				  {
					  if (face2vert[j][0] == v)
					  {
						  INDEX_3 face = new INDEX_3(face2vert[j].I1(), face2vert[j].I2(), face2vert[j].I3());
						  vert2face.Set(face, j);
					  }
					  else
					  {
						break;
					  }
				  }


				  LoopOverFaces(mesh[0], this, v, (INDEX_4 i4, int elnr, int j, bool volume, int facedir) =>
				  {
								   INDEX_3 face = new INDEX_3(i4.I1(), i4.I2(), i4.I3());
								   int facenum = vert2face.Get(face);
								   if (volume)
								   {
									   faces[elnr][j].fnr = facenum;
									   // faces[elnr][j].forient = facedir;
								   }
								   else
								   {
									   surffaces[elnr].fnr = facenum;
									   // surffaces[elnr].forient = facedir;
								   }
				  });
			  }
		});

		  /*
		  int oldnfa = face2vert.Size();
		int nfa = oldnfa;
		INDEX_3_CLOSED_HASHTABLE<int> vert2face(2*max_face_on_vertex+10);

		for (auto v : mesh.Points().Range())
		  {
		    int first_fa = nfa;

		    vert2face.DeleteData();
		  
		    for (int j = 0; j < vert2oldface[v].Size(); j++)
		      {
		        int fnr = vert2oldface[v][j];
		        INDEX_3 face (face2vert[fnr].I1(),
		                      face2vert[fnr].I2(),
		                      face2vert[fnr].I3());
		        vert2face.Set (face, fnr+1);
		      }
		  

		    for (int pass = 1; pass <= 2; pass++)
		  {

		        for (ElementIndex elnr : (*vert2element)[v])
		  {
			const Element & el = mesh[elnr];

			int nelfaces = GetNFaces (el.GetType());
			const ELEMENT_FACE * elfaces = GetFaces0 (el.GetType());

			for (int j = 0; j < nelfaces; j++)
			  if (elfaces[j][3] < 0)
		  
			{ // triangle
			  INDEX_3 face(el[elfaces[j][0]], el[elfaces[j][1]],
		                               el[elfaces[j][2]]);
		  
			  int facedir = 0;
			  if (face.I1() > face.I2())
				{ swap (face.I1(), face.I2()); facedir += 1; }
			  if (face.I2() > face.I3())
				{ swap (face.I2(), face.I3()); facedir += 2; }
			  if (face.I1() > face.I2())
				{ swap (face.I1(), face.I2()); facedir += 4; }

			  if (face.I1() != v) continue;

		                  if (pass == 1)
		                    {
		                      if (!vert2face.Used (face))
		                        {
		                          nfa++;
		                          vert2face.Set (face, nfa);
		                          INDEX_4 hface(face.I1(),face.I2(),face.I3(),0);
		                          face2vert.Append (hface);
		                        }
		                    }
		                  else
		                    {
		                      int facenum = vert2face.Get(face);
		                      faces[elnr][j].fnr = facenum-1;
		                      faces[elnr][j].forient = facedir;
		                    }
			}
		  
			  else

			{
			  // quad
			  int facenum;
			  INDEX_4Q face4(el[elfaces[j][0]], el[elfaces[j][1]],
					 el[elfaces[j][2]], el[elfaces[j][3]]);

			  int facedir = 0;
			  if (min2 (face4.I1(), face4.I2()) >
				  min2 (face4.I4(), face4.I3()))
				{  // z - flip
				  facedir += 1;
				  swap (face4.I1(), face4.I4());
				  swap (face4.I2(), face4.I3());
				}
			  if (min2 (face4.I1(), face4.I4()) >
				  min2 (face4.I2(), face4.I3()))
				{  // x - flip
				  facedir += 2;
				  swap (face4.I1(), face4.I2());
				  swap (face4.I3(), face4.I4());
				}
			  if (face4.I2() > face4.I4())
				{  // diagonal flip
				  facedir += 4;
				  swap (face4.I2(), face4.I4());
				}


			  INDEX_3 face(face4.I1(), face4.I2(), face4.I3());

			  if (face.I1() != v) continue;

			  if (vert2face.Used (face))
				{
				  facenum = vert2face.Get(face);
				}
			  else
				{
		                      if (pass == 2) cout << "hier in pass 2" << endl;
				  nfa++;
				  vert2face.Set (face, nfa);
				  facenum = nfa;

				  INDEX_4 hface(face4.I1(),face4.I2(),face4.I3(),face4.I4());
				  face2vert.Append (hface);
				}

		                  faces[elnr][j].fnr = facenum-1;
		                  faces[elnr][j].forient = facedir;
			}
		  }

		for (int j = 0; j < (*vert2surfelement)[v].Size(); j++)
		  {
			SurfaceElementIndex elnr = (*vert2surfelement)[v][j];
			const Element2d & el = mesh.SurfaceElement (elnr);

			const ELEMENT_FACE * elfaces = GetFaces1 (el.GetType());

			if (elfaces[0][3] == 0)

			  { // triangle

			int facenum;
			int facedir;

			INDEX_3 face(el.PNum(elfaces[0][0]),
					 el.PNum(elfaces[0][1]),
					 el.PNum(elfaces[0][2]));

			facedir = 0;
			if (face.I1() > face.I2())
			  {
				swap (face.I1(), face.I2());
				facedir += 1;
			  }
			if (face.I2() > face.I3())
			  {
				swap (face.I2(), face.I3());
				facedir += 2;
			  }
			if (face.I1() > face.I2())
			  {
				swap (face.I1(), face.I2());
				facedir += 4;
			  }

			if (face.I1() != v) continue;

			if (vert2face.Used (face))
			  facenum = vert2face.Get(face);
			else
			  {
				nfa++;
				vert2face.Set (face, nfa);
				facenum = nfa;

				INDEX_4 hface(face.I1(),face.I2(),face.I3(),0);
				face2vert.Append (hface);
			  }

		                surffaces[elnr].fnr = facenum-1;
		                surffaces[elnr].forient = facedir;
			  }

			else

			  {
			// quad
			int facenum;
			int facedir;

			INDEX_4Q face4(el.PNum(elfaces[0][0]),
					   el.PNum(elfaces[0][1]),
					   el.PNum(elfaces[0][2]),
					   el.PNum(elfaces[0][3]));

			facedir = 0;
			if (min2 (face4.I1(), face4.I2()) >
				min2 (face4.I4(), face4.I3()))
			  {  // z - orientation
				facedir += 1;
				swap (face4.I1(), face4.I4());
				swap (face4.I2(), face4.I3());
			  }
			if (min2 (face4.I1(), face4.I4()) >
				min2 (face4.I2(), face4.I3()))
			  {  // x - orientation
				facedir += 2;
				swap (face4.I1(), face4.I2());
				swap (face4.I3(), face4.I4());
			  }
			if (face4.I2() > face4.I4())
			  {
				facedir += 4;
				swap (face4.I2(), face4.I4());
			  }

			INDEX_3 face(face4.I1(), face4.I2(), face4.I3());
			if (face.I1() != v) continue;

			if (vert2face.Used (face))
			  facenum = vert2face.Get(face);
			else
			  {
				nfa++;
				vert2face.Set (face, nfa);
				facenum = nfa;

				INDEX_4 hface(face4.I1(),face4.I2(),face4.I3(),face4.I4());
				face2vert.Append (hface);
			  }

		                surffaces[elnr].fnr = facenum-1;
		                surffaces[elnr].forient = facedir;
			  }
		  }

		// sort faces
		if (pass == 1)
		  {
		            QuickSort (face2vert.Range(first_fa, nfa));

		            for (int j = first_fa; j < face2vert.Size(); j++)
		              {
		                if (face2vert[j][0] == v)
		                  {
		                    INDEX_3 face (face2vert[j].I1(),
		                                  face2vert[j].I2(),
		                                  face2vert[j].I3());
		                    vert2face.Set (face, j+1);
		                  }
		                else
		                  break;
		              }
		  }
		  }
	  }
		face2vert.SetAllocSize (nfa);
		  */

	// *testout << "face2vert = " << endl << face2vert << endl;

		NgProfiler.StopTimer(Update_timer2b);
		NgProfiler.StartTimer(Update_timer2c);


	face2surfel.SetSize(nfa);
	face2surfel = 0;
	for (int i = 1; i <= nse; i++)
	{
	  face2surfel.Elem(GetSurfaceElementFace(i)) = i;
	}

	/*
	  cout << "build table complete" << endl;

	  cout << "faces = " << endl;

	  cout << "face2vert = " << endl << face2vert << endl;
	  cout << "surffaces = " << endl << surffaces << endl;
	  cout << "face2surfel = " << endl << face2surfel << endl;
	*/


	surf2volelement.SetSize(nse);
	for (int i = 1; i <= nse; i++)
	{
		surf2volelement.Elem[] i = 0;
		surf2volelement.Elem[] i = 0;
	}
		tracer("Topology::Update build surf2vol", false);
	for (int i = 1; i <= ne; i++)
	{
	  for (int j = 0; j < 6; j++)
	  {
		  // int fnum = (faces.Get(i)[j]+7) / 8;
			  int fnum = faces.Get(i)[j].fnr + 1;
		  if (fnum > 0 && face2surfel.Elem(fnum))
		  {
		  int sel = face2surfel.Elem(fnum);
		  surf2volelement.Elem[] sel = surf2volelement.Elem(sel)[0];
		  surf2volelement.Elem[] sel = i;
		  }
	  }
	}
		tracer("Topology::Update build surf2vol", true);

	face2vert.SetAllocSize(face2vert.Size());

	// face table complete


#if PARALLEL
	// (*testout) << " RESET Paralleltop" << endl;
	// paralleltop.Reset ();
#endif

		tracer("Topology::Update count face_els", false);
	Array<short> face_els = new Array<short>(nfa);
	Array<short> face_surfels = new Array<short>(nfa);
	face_els = 0;
	face_surfels = 0;
		/*
	Array<int> hfaces;
	for (int i = 1; i <= ne; i++)
	  {
		GetElementFaces (i, hfaces);
		for (int j = 0; j < hfaces.Size(); j++)
		  face_els[hfaces[j]-1]++;
	  }
		*/
		netgen.GlobalMembers.ParallelForRange(tm, (uint)ne, (uint begin, uint end) =>
		{
			  Array<int> hfaces = new Array<int>();
			  for (ElementIndex ei = begin; ei < end; ei++)
			  {
				  GetElementFaces(ei + 1, hfaces);
				  foreach (var f in hfaces)
				  {
					AsAtomic(face_els[f - 1])++;
				  }
			  }
		});
	for (int i = 1; i <= nse; i++)
	{
	  face_surfels[GetSurfaceElementFace(i) - 1]++;
	}
		tracer("Topology::Update count face_els", true);


	if (ne != 0)
	{
		int cnt_err = 0;
		for (int i = 0; i < nfa; i++)
		{
		/*
		  (*testout) << "face " << i << " has " << int(face_els[i]) << " els, "
		  << int(face_surfels[i]) << " surfels, tot = "
		  << face_els[i] + face_surfels[i] << endl;
		*/
		if (face_els[i] + face_surfels[i] == 1)
		{
			cnt_err++;
#if PARALLEL
			if (ntasks > 1)
			{
			continue;
			// if ( !paralleltop.DoCoarseUpdate() ) continue;
			}
			else
#endif
			{
			(*testout) << "illegal face : " << i << "\n";
			(*testout) << "points = " << face2vert[i] << "\n";
			(*testout) << "pos = ";
			for (int j = 0; j < 4; j++)
			{
						  if (face2vert[i].I(j + 1) >= 1)
						  {
							(*testout) << mesh[(PointIndex)face2vert[i].I(j + 1)] << " ";
						  }
			}
			(*testout) << "\n";

			FlatArray<ElementIndex> vertels = GetVertexElements(face2vert[i].I(1));
			for (int k = 0; k < vertels.Size(); k++)
			{
				int[] elfaces = new int[10];
				int[] orient = new int[10];
				int nf = GetElementFaces(vertels[k] + 1, elfaces, orient);
				for (int l = 0; l < nf; l++)
				{
				  if (elfaces[l] == i)
				  {
				  // (*testout) << "is face of element " << vertels[k] << endl;

				  if (mesh.coarsemesh != null && mesh.hpelements.Size() == mesh.GetNE())
				  {
					  HPRefElement hpref_el = (mesh[0].hpelements)[mesh[vertels[k]].hp_elnr];
					  (*testout) << "coarse eleme = " << hpref_el.coarse_elnr << "\n";
				  }

				  }
				}
			}
			}
		}
		}

		if (cnt_err != 0 && ntasks == 1)
		{
		  Console.Write(cnt_err);
		  Console.Write(" elements are not matching !!!");
		  Console.Write("\n");
		}
	}
		NgProfiler.StopTimer(Update_timer2c);
	}


#if PARALLEL
	if (id != 0)
	{
	// if ( paralleltop.DoCoarseUpdate() )
	// paralleltop.UpdateCoarseGrid();
	}
#endif



	/*
	   for (i = 1; i <= ne; i++)
	   {
	   (*testout) << "Element " << i << endl;
	   (*testout) << "PNums " << endl;
	   for( int l=1;l<=8;l++) *testout << mesh.VolumeElement(i).PNum(l) << "\t";
	   *testout << endl;
	   (*testout) << "edges: " << endl;
	   for (j = 0; j < 9; j++)
	   (*testout) << edges.Elem(i)[j] << " ";
	   (*testout) << "faces: " << endl;
	   for (j = 0; j < 6; j++)m
	   (*testout) << faces.Elem(i)[j] << " ";
	   }

	   for (i = 1; i <= nse; i++)
	   {
	   (*testout) << "SElement " << i << endl;
	   (*testout) << "PNums " << endl;
	   for( int l=1;l<=4;l++) *testout << mesh.SurfaceElement(i).PNum(l) << "\t";
	   *testout << endl;
	   }
	*/
	timestamp = netgen.GlobalMembers.NextTimeStamp();
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool NeedsUpdate() const
  public bool NeedsUpdate()
  {
	  return (timestamp <= mesh.GetTimeStamp());
  }


//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetNEdges() const
  public int GetNEdges()
  {
	  return edge2vert.Size();
  }
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetNFaces() const
  public int GetNFaces()
  {
	  return face2vert.Size();
  }

  public static short GetNVertices(ELEMENT_TYPE et)
  {
	switch (et)
	{
	  case ELEMENT_TYPE.SEGMENT:
	  case ELEMENT_TYPE.SEGMENT3:
		return 2;

	  case ELEMENT_TYPE.TRIG:
	  case ELEMENT_TYPE.TRIG6:
		return 3;

	  case ELEMENT_TYPE.QUAD:
	  case ELEMENT_TYPE.QUAD6:
	  case ELEMENT_TYPE.QUAD8:
		return 4;

	  case ELEMENT_TYPE.TET:
	  case ELEMENT_TYPE.TET10:
		return 4;

	  case ELEMENT_TYPE.PYRAMID:
	  case ELEMENT_TYPE.PYRAMID13:
		return 5;

	  case ELEMENT_TYPE.PRISM:
	  case ELEMENT_TYPE.PRISM12:
	  case ELEMENT_TYPE.PRISM15:
		return 6;

	  case ELEMENT_TYPE.HEX:
	  case ELEMENT_TYPE.HEX20:
		return 8;

		// default:
		// cerr << "Ng_ME_GetNVertices, illegal element type " << et << endl;
	}
	return 0;
  }

  public static short GetNPoints(ELEMENT_TYPE et)
  {
	switch (et)
	{
	  case ELEMENT_TYPE.SEGMENT:
		return 2;
	  case ELEMENT_TYPE.SEGMENT3:
		return 3;

	  case ELEMENT_TYPE.TRIG:
		return 3;
	  case ELEMENT_TYPE.TRIG6:
		return 6;

	  case ELEMENT_TYPE.QUAD:
	  case ELEMENT_TYPE.QUAD6:
		return 4;

	  case ELEMENT_TYPE.QUAD8:
		return 8;

	  case ELEMENT_TYPE.TET:
		return 4;
	  case ELEMENT_TYPE.TET10:
		return 10;

	  case ELEMENT_TYPE.PYRAMID:
		return 5;
	  case ELEMENT_TYPE.PYRAMID13:
		return 13;

	  case ELEMENT_TYPE.PRISM:
		return 6;
	  case ELEMENT_TYPE.PRISM12:
		return 12;
	  case ELEMENT_TYPE.PRISM15:
		return 15;

	  case ELEMENT_TYPE.HEX:
		return 8;

	  case ELEMENT_TYPE.HEX20:
		return 20;
		// default:
		// cerr << "Ng_ME_GetNVertices, illegal element type " << et << endl;
	}
	return 0;
  }

  public static short GetNEdges(ELEMENT_TYPE et)
  {
	__assume(et >= ELEMENT_TYPE.SEGMENT && et <= ELEMENT_TYPE.PYRAMID13);
	switch (et)
	{
	  case ELEMENT_TYPE.SEGMENT:
	  case ELEMENT_TYPE.SEGMENT3:
		return 1;

	  case ELEMENT_TYPE.TRIG:
	  case ELEMENT_TYPE.TRIG6:
		return 3;

	  case ELEMENT_TYPE.QUAD:
	  case ELEMENT_TYPE.QUAD6:
	  case ELEMENT_TYPE.QUAD8:
		return 4;

	  case ELEMENT_TYPE.TET:
	  case ELEMENT_TYPE.TET10:
		return 6;

	  case ELEMENT_TYPE.PYRAMID:
	  case ELEMENT_TYPE.PYRAMID13:
		return 8;

	  case ELEMENT_TYPE.PRISM:
	  case ELEMENT_TYPE.PRISM12:
	  case ELEMENT_TYPE.PRISM15:
		return 9;

	  case ELEMENT_TYPE.HEX:
	  case ELEMENT_TYPE.HEX20:
		return 12;
	  default:
		return 0;
		// default:
		// cerr << "Ng_ME_GetNEdges, illegal element type " << et << endl;
	}
	// return 0;
  }

  public static short GetNFaces(ELEMENT_TYPE et)
  {
	__assume(et >= ELEMENT_TYPE.SEGMENT && et <= ELEMENT_TYPE.PYRAMID13);
	switch (et)
	{
	  case ELEMENT_TYPE.SEGMENT:
	  case ELEMENT_TYPE.SEGMENT3:
		return 0;

	  case ELEMENT_TYPE.TRIG:
	  case ELEMENT_TYPE.TRIG6:
		return 1;

	  case ELEMENT_TYPE.QUAD:
	  case ELEMENT_TYPE.QUAD6:
	  case ELEMENT_TYPE.QUAD8:
		return 1;

	  case ELEMENT_TYPE.TET:
	  case ELEMENT_TYPE.TET10:
		return 4;

	  case ELEMENT_TYPE.PYRAMID:
	  case ELEMENT_TYPE.PYRAMID13:
		return 5;

	  case ELEMENT_TYPE.PRISM:
	  case ELEMENT_TYPE.PRISM12:
	  case ELEMENT_TYPE.PRISM15:
		return 5;

	  case ELEMENT_TYPE.HEX:
	  case ELEMENT_TYPE.HEX20:
		return 6;

	  default:
		return -99;
		// default:
		// cerr << "Ng_ME_GetNVertices, illegal element type " << et << endl;
	}
  }

//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private static Point3d[] GetVertices_segm_points = {Point3d(1, 0, 0), Point3d(0, 0, 0)};
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private static Point3d[] GetVertices_trig_points = {Point3d(1, 0, 0), Point3d(0, 1, 0), Point3d(0, 0, 0)};
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private static Point3d[] GetVertices_quad_points = {Point3d(0, 0, 0), Point3d(1, 0, 0), Point3d(1, 1, 0), Point3d(0, 1, 0)};
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private static Point3d[] GetVertices_tet_points = {Point3d(1, 0, 0), Point3d(0, 1, 0), Point3d(0, 0, 1), Point3d(0, 0, 0)};
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private static Point3d[] GetVertices_pyramid_points = {Point3d(0, 0, 0), Point3d(1, 0, 0), Point3d(1, 1, 0), Point3d(0, 1, 0), Point3d(0, 0, 1 - 1e-7)};
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private static Point3d[] GetVertices_prism_points = {Point3d(1, 0, 0), Point3d(0, 1, 0), Point3d(0, 0, 0), Point3d(1, 0, 1), Point3d(0, 1, 1), Point3d(0, 0, 1)};
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private static Point3d[] GetVertices_hex_points = {Point3d(0, 0, 0), Point3d(1, 0, 0), Point3d(1, 1, 0), Point3d(0, 1, 0), Point3d(0, 0, 1), Point3d(1, 0, 1), Point3d(1, 1, 1), Point3d(0, 1, 1)};

  public static Point3d GetVertices(ELEMENT_TYPE et)
  {
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static Point3d segm_points [] = { Point3d(1, 0, 0), Point3d(0, 0, 0) };

//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static Point3d trig_points [] = { Point3d(1, 0, 0), Point3d(0, 1, 0), Point3d(0, 0, 0) };

//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static Point3d quad_points [] = { Point3d(0, 0, 0), Point3d(1, 0, 0), Point3d(1, 1, 0), Point3d(0, 1, 0) };

//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static Point3d tet_points [] = { Point3d(1, 0, 0), Point3d(0, 1, 0), Point3d(0, 0, 1), Point3d(0, 0, 0) };

//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static Point3d pyramid_points [] = { Point3d(0, 0, 0), Point3d(1, 0, 0), Point3d(1, 1, 0), Point3d(0, 1, 0), Point3d(0, 0, 1-1e-7)};

//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static Point3d prism_points[] = { Point3d(1, 0, 0), Point3d(0, 1, 0), Point3d(0, 0, 0), Point3d(1, 0, 1), Point3d(0, 1, 1), Point3d(0, 0, 1) };


//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static Point3d hex_points [] = { Point3d(0, 0, 0), Point3d(1, 0, 0), Point3d(1, 1, 0), Point3d(0, 1, 0), Point3d(0, 0, 1), Point3d(1, 0, 1), Point3d(1, 1, 1), Point3d(0, 1, 1) };


	switch (et)
	{
	  case ELEMENT_TYPE.SEGMENT:
	  case ELEMENT_TYPE.SEGMENT3:
	return GetVertices_segm_points;

	  case ELEMENT_TYPE.TRIG:
	  case ELEMENT_TYPE.TRIG6:
	return GetVertices_trig_points;

	  case ELEMENT_TYPE.QUAD:
	  case ELEMENT_TYPE.QUAD6:
	  case ELEMENT_TYPE.QUAD8:
	return GetVertices_quad_points;

	  case ELEMENT_TYPE.TET:
	  case ELEMENT_TYPE.TET10:
	return GetVertices_tet_points;

	  case ELEMENT_TYPE.PYRAMID:
	return GetVertices_pyramid_points;

	  case ELEMENT_TYPE.PRISM:
	  case ELEMENT_TYPE.PRISM12:
	return GetVertices_prism_points;

	  case ELEMENT_TYPE.HEX:
	return GetVertices_hex_points;
	  default:
	cerr << "Ng_ME_GetVertices, illegal element type " << et << "\n";
	break;
	}
	return null;
  }

//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private static ELEMENT_EDGE[] GetEdges1_segm_edges =
  {
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {1, 2}
	  }
  };
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private static ELEMENT_EDGE[] GetEdges1_trig_edges =
  {
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {3, 1}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {2, 3}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {1, 2}
	  }
  };
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private static ELEMENT_EDGE[] GetEdges1_quad_edges =
  {
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {1, 2}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {3, 4}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {4, 1}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {2, 3}
	  }
  };
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private static ELEMENT_EDGE[] GetEdges1_tet_edges =
  {
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {4, 1}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {4, 2}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {4, 3}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {1, 2}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {1, 3}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {2, 3}
	  }
  };
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private static ELEMENT_EDGE[] GetEdges1_prism_edges =
  {
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {3, 1}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {1, 2}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {3, 2}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {6, 4}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {4, 5}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {6, 5}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {3, 6}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {1, 4}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {2, 5}
	  }
  };
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private static ELEMENT_EDGE[] GetEdges1_pyramid_edges =
  {
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {1, 2}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {2, 3}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {1, 4}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {4, 3}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {1, 5}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {2, 5}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {3, 5}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {4, 5}
	  }
  };
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private static ELEMENT_EDGE[] GetEdges1_hex_edges =
  {
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {1, 2}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {3, 4}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {4, 1}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {2, 3}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {5, 6}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {7, 8}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {8, 5}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {6, 7}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {1, 5}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {2, 6}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {3, 7}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {4, 8}
	  }
  };

  public static ELEMENT_EDGE GetEdges1(ELEMENT_TYPE et)
  {
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static ELEMENT_EDGE segm_edges[1] = { { 1, 2 }};

//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static ELEMENT_EDGE trig_edges[3] = { { 3, 1 }, { 2, 3 }, { 1, 2 }};

//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static ELEMENT_EDGE quad_edges[4] = { { 1, 2 }, { 3, 4 }, { 4, 1 }, { 2, 3 }};


//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static ELEMENT_EDGE tet_edges[6] = { { 4, 1 }, { 4, 2 }, { 4, 3 }, { 1, 2 }, { 1, 3 }, { 2, 3 }};

//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static ELEMENT_EDGE prism_edges[9] = { { 3, 1 }, { 1, 2 }, { 3, 2 }, { 6, 4 }, { 4, 5 }, { 6, 5 }, { 3, 6 }, { 1, 4 }, { 2, 5 }};

//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static ELEMENT_EDGE pyramid_edges[8] = { { 1, 2 }, { 2, 3 }, { 1, 4 }, { 4, 3 }, { 1, 5 }, { 2, 5 }, { 3, 5 }, { 4, 5 }};

//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static ELEMENT_EDGE hex_edges[12] = { { 1, 2 }, { 3, 4 }, { 4, 1 }, { 2, 3 }, { 5, 6 }, { 7, 8 }, { 8, 5 }, { 6, 7 }, { 1, 5 }, { 2, 6 }, { 3, 7 }, { 4, 8 }};

	switch (et)
	{
	  case ELEMENT_TYPE.SEGMENT:
	  case ELEMENT_TYPE.SEGMENT3:
		return GetEdges1_segm_edges;

	  case ELEMENT_TYPE.TRIG:
	  case ELEMENT_TYPE.TRIG6:
		return GetEdges1_trig_edges;

	  case ELEMENT_TYPE.QUAD:
	  case ELEMENT_TYPE.QUAD6:
	  case ELEMENT_TYPE.QUAD8:
		return GetEdges1_quad_edges;

	  case ELEMENT_TYPE.TET:
	  case ELEMENT_TYPE.TET10:
		return GetEdges1_tet_edges;

	  case ELEMENT_TYPE.PYRAMID:
	  case ELEMENT_TYPE.PYRAMID13:
		return GetEdges1_pyramid_edges;

	  case ELEMENT_TYPE.PRISM:
	  case ELEMENT_TYPE.PRISM12:
	  case ELEMENT_TYPE.PRISM15:
		return GetEdges1_prism_edges;

	  case ELEMENT_TYPE.HEX:
	  case ELEMENT_TYPE.HEX20:
		return GetEdges1_hex_edges;
		// default:
		// cerr << "Ng_ME_GetEdges, illegal element type " << et << endl;
	}
	 return null;
  }

//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private static ELEMENT_EDGE[] GetEdges0_segm_edges =
  {
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {0, 1}
	  }
  };
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private static ELEMENT_EDGE[] GetEdges0_trig_edges =
  {
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {2, 0}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {1, 2}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {0, 1}
	  }
  };
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private static ELEMENT_EDGE[] GetEdges0_quad_edges =
  {
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {0, 1}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {2, 3}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {3, 0}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {1, 2}
	  }
  };
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private static ELEMENT_EDGE[] GetEdges0_tet_edges =
  {
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {3, 0}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {3, 1}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {3, 2}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {0, 1}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {0, 2}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {1, 2}
	  }
  };
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private static ELEMENT_EDGE[] GetEdges0_prism_edges =
  {
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {2, 0}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {0, 1}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {2, 1}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {5, 3}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {3, 4}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {5, 4}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {2, 5}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {0, 3}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {1, 4}
	  }
  };
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private static ELEMENT_EDGE[] GetEdges0_pyramid_edges =
  {
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {0, 1}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {1, 2}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {0, 3}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {3, 2}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {0, 4}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {1, 4}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {2, 4}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {3, 4}
	  }
  };
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private static ELEMENT_EDGE[] GetEdges0_hex_edges =
  {
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {0, 1}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {2, 3}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {3, 0}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {1, 2}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {4, 5}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {6, 7}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {7, 4}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {5, 6}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {0, 4}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {1, 5}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {2, 6}
	  },
	  new ELEMENT_EDGE()
	  {
		  vals = new[] {3, 7}
	  }
  };

  public static ELEMENT_EDGE GetEdges0(ELEMENT_TYPE et)
  {
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static ELEMENT_EDGE segm_edges[1] = { { 0, 1 }};

//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static ELEMENT_EDGE trig_edges[3] = { { 2, 0 }, { 1, 2 }, { 0, 1 }};

//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static ELEMENT_EDGE quad_edges[4] = { { 0, 1 }, { 2, 3 }, { 3, 0 }, { 1, 2 }};


//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static ELEMENT_EDGE tet_edges[6] = { { 3, 0 }, { 3, 1 }, { 3, 2 }, { 0, 1 }, { 0, 2 }, { 1, 2 }};

//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static ELEMENT_EDGE prism_edges[9] = { { 2, 0 }, { 0, 1 }, { 2, 1 }, { 5, 3 }, { 3, 4 }, { 5, 4 }, { 2, 5 }, { 0, 3 }, { 1, 4 }};

//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static ELEMENT_EDGE pyramid_edges[8] = { { 0, 1 }, { 1, 2 }, { 0, 3 }, { 3, 2 }, { 0, 4 }, { 1, 4 }, { 2, 4 }, { 3, 4 }};

//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static ELEMENT_EDGE hex_edges[12] = { { 0, 1 }, { 2, 3 }, { 3, 0 }, { 1, 2 }, { 4, 5 }, { 6, 7 }, { 7, 4 }, { 5, 6 }, { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 }};

	switch (et)
	{
	  case ELEMENT_TYPE.SEGMENT:
	  case ELEMENT_TYPE.SEGMENT3:
		return GetEdges0_segm_edges;

	  case ELEMENT_TYPE.TRIG:
	  case ELEMENT_TYPE.TRIG6:
		return GetEdges0_trig_edges;

	  case ELEMENT_TYPE.QUAD:
	  case ELEMENT_TYPE.QUAD6:
	  case ELEMENT_TYPE.QUAD8:
		return GetEdges0_quad_edges;

	  case ELEMENT_TYPE.TET:
	  case ELEMENT_TYPE.TET10:
		return GetEdges0_tet_edges;

	  case ELEMENT_TYPE.PYRAMID:
	  case ELEMENT_TYPE.PYRAMID13:
		return GetEdges0_pyramid_edges;

	  case ELEMENT_TYPE.PRISM:
	  case ELEMENT_TYPE.PRISM12:
	  case ELEMENT_TYPE.PRISM15:
		return GetEdges0_prism_edges;

	  case ELEMENT_TYPE.HEX:
	  case ELEMENT_TYPE.HEX20:
		return GetEdges0_hex_edges;
		// default:
		// cerr << "Ng_ME_GetEdges, illegal element type " << et << endl;
	}
	 return null;
  }

  public static ELEMENT_FACE[] GetFaces1(ELEMENT_TYPE et)
  {
	ELEMENT_FACE[] trig_faces =
	{
		new ELEMENT_FACE()
		{
			vals = new[] {1, 2, 3, 0}
		}
	};
	ELEMENT_FACE[] quad_faces =
	{
		new ELEMENT_FACE()
		{
			vals = new[] {1, 2, 3, 4}
		}
	};

	ELEMENT_FACE[] tet_faces =
	{
		new ELEMENT_FACE()
		{
			vals = new[] {4, 2, 3, 0}
		},
		new ELEMENT_FACE()
		{
			vals = new[] {4, 3, 1, 0}
		},
		new ELEMENT_FACE()
		{
			vals = new[] {4, 1, 2, 0}
		},
		new ELEMENT_FACE()
		{
			vals = new[] {1, 3, 2, 0}
		}
	};

	ELEMENT_FACE[] prism_faces =
	{
		new ELEMENT_FACE()
		{
			vals = new[] {1, 3, 2, 0}
		},
		new ELEMENT_FACE()
		{
			vals = new[] {4, 5, 6, 0}
		},
		new ELEMENT_FACE()
		{
			vals = new[] {3, 1, 4, 6}
		},
		new ELEMENT_FACE()
		{
			vals = new[] {1, 2, 5, 4}
		},
		new ELEMENT_FACE()
		{
			vals = new[] {2, 3, 6, 5}
		}
	};

	ELEMENT_FACE[] pyramid_faces =
	{
		new ELEMENT_FACE()
		{
			vals = new[] {1, 2, 5, 0}
		},
		new ELEMENT_FACE()
		{
			vals = new[] {2, 3, 5, 0}
		},
		new ELEMENT_FACE()
		{
			vals = new[] {3, 4, 5, 0}
		},
		new ELEMENT_FACE()
		{
			vals = new[] {4, 1, 5, 0}
		},
		new ELEMENT_FACE()
		{
			vals = new[] {1, 4, 3, 2}
		}
	};

	ELEMENT_FACE[] hex_faces =
	{
		new ELEMENT_FACE()
		{
			vals = new[] {1, 4, 3, 2}
		},
		new ELEMENT_FACE()
		{
			vals = new[] {5, 6, 7, 8}
		},
		new ELEMENT_FACE()
		{
			vals = new[] {1, 2, 6, 5}
		},
		new ELEMENT_FACE()
		{
			vals = new[] {2, 3, 7, 6}
		},
		new ELEMENT_FACE()
		{
			vals = new[] {3, 4, 8, 7}
		},
		new ELEMENT_FACE()
		{
			vals = new[] {4, 1, 5, 8}
		}
	};



	switch (et)
	{
	  case ELEMENT_TYPE.TRIG:
	  case ELEMENT_TYPE.TRIG6:
		return trig_faces;

	  case ELEMENT_TYPE.QUAD:
	  case ELEMENT_TYPE.QUAD6:
	  case ELEMENT_TYPE.QUAD8:
		return quad_faces;


	  case ELEMENT_TYPE.TET:
	  case ELEMENT_TYPE.TET10:
		return tet_faces;

	  case ELEMENT_TYPE.PRISM:
	  case ELEMENT_TYPE.PRISM12:
	  case ELEMENT_TYPE.PRISM15:
		return prism_faces;

	  case ELEMENT_TYPE.PYRAMID:
	  case ELEMENT_TYPE.PYRAMID13:
		return pyramid_faces;

	  case ELEMENT_TYPE.SEGMENT:
	  case ELEMENT_TYPE.SEGMENT3:

	  case ELEMENT_TYPE.HEX:
	  case ELEMENT_TYPE.HEX20:
		return hex_faces;

		// default:
		// cerr << "Ng_ME_GetVertices, illegal element type " << et << endl;
	}
	return null;
  }

  public static ELEMENT_FACE[] GetFaces0(ELEMENT_TYPE et)
  {
	ELEMENT_FACE[] trig_faces =
	{
		new ELEMENT_FACE()
		{
			vals = new[] {0, 1, 2, -1}
		}
	};
	ELEMENT_FACE[] quad_faces =
	{
		new ELEMENT_FACE()
		{
			vals = new[] {0, 1, 2, 3}
		}
	};

	ELEMENT_FACE[] tet_faces =
	{
		new ELEMENT_FACE()
		{
			vals = new[] {3, 1, 2, -1}
		},
		new ELEMENT_FACE()
		{
			vals = new[] {3, 2, 0, -1}
		},
		new ELEMENT_FACE()
		{
			vals = new[] {3, 0, 1, -1}
		},
		new ELEMENT_FACE()
		{
			vals = new[] {0, 2, 1, -1}
		}
	};

	ELEMENT_FACE[] prism_faces =
	{
		new ELEMENT_FACE()
		{
			vals = new[] {0, 2, 1, -1}
		},
		new ELEMENT_FACE()
		{
			vals = new[] {3, 4, 5, -1}
		},
		new ELEMENT_FACE()
		{
			vals = new[] {2, 0, 3, 5}
		},
		new ELEMENT_FACE()
		{
			vals = new[] {0, 1, 4, 3}
		},
		new ELEMENT_FACE()
		{
			vals = new[] {1, 2, 5, 4}
		}
	};

	ELEMENT_FACE[] pyramid_faces =
	{
		new ELEMENT_FACE()
		{
			vals = new[] {0, 1, 4, -1}
		},
		new ELEMENT_FACE()
		{
			vals = new[] {1, 2, 4, -1}
		},
		new ELEMENT_FACE()
		{
			vals = new[] {2, 3, 4, -1}
		},
		new ELEMENT_FACE()
		{
			vals = new[] {3, 0, 4, -1}
		},
		new ELEMENT_FACE()
		{
			vals = new[] {0, 3, 2, 1}
		}
	};

	ELEMENT_FACE[] hex_faces =
	{
		new ELEMENT_FACE()
		{
			vals = new[] {0, 3, 2, 1}
		},
		new ELEMENT_FACE()
		{
			vals = new[] {4, 5, 6, 7}
		},
		new ELEMENT_FACE()
		{
			vals = new[] {0, 1, 5, 4}
		},
		new ELEMENT_FACE()
		{
			vals = new[] {1, 2, 6, 5}
		},
		new ELEMENT_FACE()
		{
			vals = new[] {2, 3, 7, 6}
		},
		new ELEMENT_FACE()
		{
			vals = new[] {3, 0, 4, 7}
		}
	};



	switch (et)
	{
	  case ELEMENT_TYPE.TRIG:
	  case ELEMENT_TYPE.TRIG6:
		return trig_faces;

	  case ELEMENT_TYPE.QUAD:
	  case ELEMENT_TYPE.QUAD6:
	  case ELEMENT_TYPE.QUAD8:
		return quad_faces;


	  case ELEMENT_TYPE.TET:
	  case ELEMENT_TYPE.TET10:
		return tet_faces;

	  case ELEMENT_TYPE.PRISM:
	  case ELEMENT_TYPE.PRISM12:
	  case ELEMENT_TYPE.PRISM15:
		return prism_faces;

	  case ELEMENT_TYPE.PYRAMID:
	  case ELEMENT_TYPE.PYRAMID13:
		return pyramid_faces;

	  case ELEMENT_TYPE.SEGMENT:
	  case ELEMENT_TYPE.SEGMENT3:

	  case ELEMENT_TYPE.HEX:
	  case ELEMENT_TYPE.HEX20:
		return hex_faces;

		// default:
		// cerr << "Ng_ME_GetVertices, illegal element type " << et << endl;
	}
	return null;
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetSegmentEdge(int segnr) const
  public int GetSegmentEdge(int segnr)
  {
	  return segedges[segnr - 1].nr + 1;
  }
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetEdge(SegmentIndex segnr) const
  public int GetEdge(SegmentIndex segnr)
  {
	  return segedges[segnr].nr;
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetSegmentEdge(int segnr, int & enr, int & orient) const
  public void GetSegmentEdge(int segnr, ref int enr, ref int orient)
  {
	enr = segedges.Get(segnr).nr + 1;
	// orient = segedges.Get(segnr).orient;
	orient = GetSegmentEdgeOrientation(segnr);
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetElementEdges(int elnr, Array<int> & eledges) const
  public void GetElementEdges(int elnr, Array<int> eledges)
  {
	int ned = GetNEdges(mesh.VolumeElement(elnr).GetType());
	eledges.SetSize(ned);
	for (int i = 0; i < ned; i++)
	{
	  eledges[i] = edges.Get(elnr)[i].nr + 1;
	}
	  // eledges[i] = abs (edges.Get(elnr)[i]);
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetElementFaces(int elnr, Array<int> & elfaces, bool withorientation = false) const
  public void GetElementFaces(int elnr, Array<int> elfaces, bool withorientation = false)
  {
	int nfa = GetNFaces(mesh.VolumeElement(elnr).GetType());
	elfaces.SetSize(nfa);

	if (!withorientation)
	{

	  for (int i = 1; i <= nfa; i++)
	  {
	  // elfaces.Elem(i) = (faces.Get(elnr)[i-1]-1) / 8 + 1;
		  elfaces.Elem(i) = faces.Get(elnr)[i - 1].fnr + 1;
	  }
	}

	else
	{
		cerr << "GetElementFaces with orientation currently not supported" << "\n";
		/*
		  for (int i = 1; i <= nfa; i++)
		  {
	  elfaces.Elem(i) = (faces.Get(elnr)[i-1]-1) / 8 + 1;
	  int orient = (faces.Get(elnr)[i-1]-1) % 8;
	  if(orient == 1 || orient == 2 || orient == 4 || orient == 7)
		  elfaces.Elem(i) *= -1;
		  }
		*/
	}
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetElementEdgeOrientations(int elnr, Array<int> & eorient) const
  public void GetElementEdgeOrientations(int elnr, Array<int> eorient)
  {
	int ned = GetNEdges(mesh.VolumeElement(elnr).GetType());
	eorient.SetSize(ned);
	for (int i = 1; i <= ned; i++)
	{
	  // eorient.Elem(i) = (edges.Get(elnr)[i-1] > 0) ? 1 : -1;
	  // eorient.Elem(i) = (edges.Get(elnr)[i-1].orient) ? -1 : 1;
	  eorient.Elem(i) = GetElementEdgeOrientation(elnr, i - 1) != 0 ? -1 : 1;
	}
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetElementFaceOrientations(int elnr, Array<int> & forient) const
  public void GetElementFaceOrientations(int elnr, Array<int> forient)
  {
	int nfa = GetNFaces(mesh.VolumeElement(elnr).GetType());
	forient.SetSize(nfa);
	for (int i = 1; i <= nfa; i++)
	{
	  // forient.Elem(i) = faces.Get(elnr)[i-1].forient;
	  // forient.Elem(i) = (faces.Get(elnr)[i-1]-1) % 8;
	  forient.Elem(i) = GetElementFaceOrientation(elnr, i - 1);
	}
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetElementEdges(int elnr, int * eledges, int * orient) const
  public int GetElementEdges(int elnr, int[] eledges, int[] orient)
  {
	//  int ned = GetNEdges (mesh.VolumeElement(elnr).GetType());

	if (mesh.GetDimension() == 3 || true)
	{
		if (orient != 0)
		{
		for (int i = 0; i < 12; i++)
		{
				/*
		if (!edges.Get(elnr)[i]) return i;
		eledges[i] = abs (edges.Get(elnr)[i]);
		orient[i] = (edges.Get(elnr)[i] > 0 ) ? 1 : -1;
				*/
				if (edges.Get(elnr)[i].nr == -1)
				{
					return i;
				}
				eledges[i] = edges.Get(elnr)[i].nr + 1;
		// orient[i] = edges.Get(elnr)[i].orient ? -1 : 1;
				orient[i] = GetElementEdgeOrientation(elnr, i) != 0 ? -1 : 1;
		}
		}
	else
	{
		for (int i = 0; i < 12; i++)
		{
		// if (!edges.Get(elnr)[i]) return i;
		// eledges[i] = abs (edges.Get(elnr)[i]);
				if (edges.Get(elnr)[i].nr == -1)
				{
					return i;
				}
				eledges[i] = edges.Get(elnr)[i].nr + 1;

		}
	}
	return 12;
	}
	else
	{
	throw new Exception("rethink implementation");
	/*
	  if (orient)
	  {
	  for (i = 0; i < 4; i++)
	  {
	  if (!surfedges.Get(elnr)[i]) return i;
	  eledges[i] = abs (surfedges.Get(elnr)[i]);
	  orient[i] = (surfedges.Get(elnr)[i] > 0 ) ? 1 : -1;
	  }
	  }
	  else
	  {
	  if (!surfedges.Get(elnr)[i]) return i;
	  for (i = 0; i < 4; i++)
	  eledges[i] = abs (surfedges.Get(elnr)[i]);
	  }
	*/
	return 4;
	//      return GetSurfaceElementEdges (elnr, eledges, orient);
	}
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetElementFaces(int elnr, int * elfaces, int * orient) const
  public int GetElementFaces(int elnr, int[] elfaces, int[] orient)
  {
	//  int nfa = GetNFaces (mesh.VolumeElement(elnr).GetType());
	if (orient != 0)
	{
	for (int i = 0; i < 6; i++)
	{
			/*
		if (!faces.Get(elnr)[i]) return i;
		elfaces[i] = (faces.Get(elnr)[i]-1) / 8 + 1;
		orient[i] = (faces.Get(elnr)[i]-1) % 8;
			*/
		if (faces.Get(elnr)[i].fnr == -1)
		{
			return i;
		}
		elfaces[i] = faces.Get(elnr)[i].fnr + 1;
		// orient[i] = faces.Get(elnr)[i].forient;
			orient[i] = GetElementFaceOrientation(elnr, i);
	}
	}
	else
	{
	for (int i = 0; i < 6; i++)
	{
		// if (!faces.Get(elnr)[i]) return i;
		// elfaces[i] = (faces.Get(elnr)[i]-1) / 8 + 1;
		if (faces.Get(elnr)[i].fnr == -1)
		{
			return i;
		}
		elfaces[i] = faces.Get(elnr)[i].fnr + 1;
	}
	}
	return 6;
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetElementEdgeOrientation(int elnr, int locedgenr) const
  public int GetElementEdgeOrientation(int elnr, int locedgenr)
  {
	Element el = mesh.VolumeElement(elnr);
	ELEMENT_EDGE[] eledges = MeshTopology.GetEdges0(el.GetType());

	int k = locedgenr;
	INDEX_2 edge = new INDEX_2(el[eledges[k][0]], el[eledges[k][1]]);
	int edgedir = (edge.I1() > edge.I2());
	return edgedir;
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetElementFaceOrientation(int elnr, int locfacenr) const
  public int GetElementFaceOrientation(int elnr, int locfacenr)
  {
	Element el = mesh.VolumeElement(elnr);

	ELEMENT_FACE[] elfaces = MeshTopology.GetFaces0(el.GetType());

	int j = locfacenr;
	if (elfaces[j][3] < 0)
	{ // triangle
		INDEX_4 face = new INDEX_4(el[elfaces[j][0]], el[elfaces[j][1]], el[elfaces[j][2]], 0);

		int facedir = 0;
		if (face.I1() > face.I2())
		{
			  swap(face.I1(), face.I2());
			  facedir += 1;
		}
		if (face.I2() > face.I3())
		{
			  swap(face.I2(), face.I3());
			  facedir += 2;
		}
		if (face.I1() > face.I2())
		{
			  swap(face.I1(), face.I2());
			  facedir += 4;
		}

		return facedir;
	}
	else
	{
		// quad
		// int facenum;
		INDEX_4 face4 = new INDEX_4(el[elfaces[j][0]], el[elfaces[j][1]], el[elfaces[j][2]], el[elfaces[j][3]]);

		int facedir = 0;
		if (netgen.GlobalMembers.min2(face4.I1(), face4.I2()) > netgen.GlobalMembers.min2(face4.I4(), face4.I3()))
		{ // z - flip
			facedir += 1;
			swap(face4.I1(), face4.I4());
			swap(face4.I2(), face4.I3());
		}
		if (netgen.GlobalMembers.min2(face4.I1(), face4.I4()) > netgen.GlobalMembers.min2(face4.I2(), face4.I3()))
		{ // x - flip
			facedir += 2;
			swap(face4.I1(), face4.I2());
			swap(face4.I3(), face4.I4());
		}
		if (face4.I2() > face4.I4())
		{ // diagonal flip
			facedir += 4;
			swap(face4.I2(), face4.I4());
		}

		return facedir;
	}
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetSurfaceElementEdgeOrientation(int elnr, int locedgenr) const
  public int GetSurfaceElementEdgeOrientation(int elnr, int locedgenr)
  {
	Element2d el = mesh.SurfaceElement(elnr);
	ELEMENT_EDGE[] eledges = MeshTopology.GetEdges0(el.GetType());

	int k = locedgenr;
	INDEX_2 edge = new INDEX_2(el[eledges[k][0]], el[eledges[k][1]]);
	int edgedir = (edge.I1() > edge.I2());
	return edgedir;
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetSurfaceElementFaceOrientation2(int elnr) const
  public int GetSurfaceElementFaceOrientation2(int elnr)
  {
	Element2d el = mesh.SurfaceElement(elnr);

	ELEMENT_FACE[] elfaces = MeshTopology.GetFaces0(el.GetType());

	int j = 0;
	if (elfaces[j][3] < 0)
	{ // triangle
		INDEX_4 face = new INDEX_4(el[elfaces[j][0]], el[elfaces[j][1]], el[elfaces[j][2]], 0);

		int facedir = 0;
		if (face.I1() > face.I2())
		{
			  swap(face.I1(), face.I2());
			  facedir += 1;
		}
		if (face.I2() > face.I3())
		{
			  swap(face.I2(), face.I3());
			  facedir += 2;
		}
		if (face.I1() > face.I2())
		{
			  swap(face.I1(), face.I2());
			  facedir += 4;
		}

		return facedir;
	}
	else
	{
		// quad
		// int facenum;
		INDEX_4 face4 = new INDEX_4(el[elfaces[j][0]], el[elfaces[j][1]], el[elfaces[j][2]], el[elfaces[j][3]]);

		int facedir = 0;
		if (netgen.GlobalMembers.min2(face4.I1(), face4.I2()) > netgen.GlobalMembers.min2(face4.I4(), face4.I3()))
		{ // z - flip
			facedir += 1;
			swap(face4.I1(), face4.I4());
			swap(face4.I2(), face4.I3());
		}
		if (netgen.GlobalMembers.min2(face4.I1(), face4.I4()) > netgen.GlobalMembers.min2(face4.I2(), face4.I3()))
		{ // x - flip
			facedir += 2;
			swap(face4.I1(), face4.I2());
			swap(face4.I3(), face4.I4());
		}
		if (face4.I2() > face4.I4())
		{ // diagonal flip
			facedir += 4;
			swap(face4.I2(), face4.I4());
		}

		return facedir;
	}
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetSegmentEdgeOrientation(int elnr) const
  public int GetSegmentEdgeOrientation(int elnr)
  {
	Segment el = mesh.LineSegment(elnr);
	ELEMENT_EDGE[] eledges = MeshTopology.GetEdges0(el.GetType());

	int k = 0;
	INDEX_2 edge = new INDEX_2(el[eledges[k][0]], el[eledges[k][1]]);
	int edgedir = (edge.I1() > edge.I2());
	return edgedir;
  }


//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetFaceVertices(int fnr, Array<int> & vertices) const
  public void GetFaceVertices(int fnr, Array<int> vertices)
  {
	vertices.SetSize(4);
	for (int i = 0; i < 4; i++)
	{
	  vertices[i] = face2vert.Get(fnr)[i];
	}
	if (vertices[3] == 0)
	{
	  vertices.SetSize(3);
	}
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetFaceVertices(int fnr, int * vertices) const
  public void GetFaceVertices(int fnr, int[] vertices)
  {
	for (int i = 0; i <= 3; i++)
	{
	  vertices[i] = face2vert.Get(fnr)[i];
	}
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetEdgeVertices(int ednr, int & v1, int & v2) const
  public void GetEdgeVertices(int ednr, ref int v1, ref int v2)
  {
	// cout << "id = " << id << "getedgevertices, ednr = " << ednr << ", ned = " << edge2vert.Size() << "&v1 = " << &v1 << endl;
	if (ednr < 1 || ednr > edge2vert.Size())
	{
	  cerr << "illegal edge nr: " << ednr << ", numedges = " << edge2vert.Size() << " id = " << id << "\n";
	}
	v1 = edge2vert.Get(ednr)[0];
	v2 = edge2vert.Get(ednr)[1];
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetEdgeVertices(int ednr, PointIndex & v1, PointIndex & v2) const
  public void GetEdgeVertices(int ednr, ref PointIndex v1, ref PointIndex v2)
  {
	v1 = edge2vert.Get(ednr)[0];
	v2 = edge2vert.Get(ednr)[1];
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const int * GetEdgeVerticesPtr(int enr) const
//C++ TO C# CONVERTER WARNING: C# has no equivalent to methods returning pointers to value types:
  public int GetEdgeVerticesPtr(int enr)
  {
	  return edge2vert[enr][0];
  }
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const int * GetFaceVerticesPtr(int fnr) const
//C++ TO C# CONVERTER WARNING: C# has no equivalent to methods returning pointers to value types:
  public int GetFaceVerticesPtr(int fnr)
  {
	  return face2vert[fnr][0];
  }
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetFaceEdges(int fnr, Array<int> & fedges, bool withorientation = false) const
  public void GetFaceEdges(int fnr, Array<int> fedges, bool withorientation = false)
  {
	ArrayMem<int,4> pi = new ArrayMem<int,4>(4);
	ArrayMem<int,12> eledges = new ArrayMem<int,12>();

	fedges.SetSize(0);
	GetFaceVertices(fnr, pi);

	// Sort Edges according to global vertex numbers
	// e1 = fmax, f2
	// e2 = fmax, f1
	// e3 = op e1(f2,f3)
	// e4 = op e2(f1,f3)

	/*  ArrayMem<int,4> fp;
	fp[0] = pi[0];
	for(int k=1;k<pi.Size();k++)
	if(fp[k]>fp[0]) swap(fp[k],fp[0]);

	fp[1] = fp[0]+ */


	//  GetVertexElements (pi[0], els);
	FlatArray<ElementIndex> els = GetVertexElements(pi[0]);

	// find one element having all vertices of the face
	for (int i = 0; i < els.Size(); i++)
	{
	Element el = mesh[els[i]];
	int nref_faces = GetNFaces(el.GetType());
	ELEMENT_FACE[] ref_faces = GetFaces1(el.GetType());
	int nfa_ref_edges = GetNEdges(GetFaceType(fnr));

	int cntv = 0;
	int fa = -1;
	for (int m = 0;m < nref_faces;m++)
	{
		cntv = 0;
		for (int j = 0;j<nfa_ref_edges && ref_faces[m][j]>0;j++)
		{
		  for (int k = 0;k < pi.Size();k++)
		  {
		  if (el[ref_faces[m][j] - 1] == pi[k])
		  {
			cntv++;
		  }
		  }
		}
		if (cntv == pi.Size())
		{
		fa = m;
		break;
		}
	}

	if (fa >= 0)
	{
		ELEMENT_EDGE[] fa_ref_edges = GetEdges1(GetFaceType(fnr));
		fedges.SetSize(nfa_ref_edges);
		GetElementEdges(els[i] + 1, eledges);

		for (int j = 0; j < eledges.Size(); j++)
		{
		int vi1;
		int vi2;
		GetEdgeVertices(eledges[j], vi1, vi2);

		bool has1 = false;
		bool has2 = false;
		for (int k = 0; k < pi.Size(); k++)
		{
			if (vi1 == pi[k])
			{
				has1 = true;
			}
			if (vi2 == pi[k])
			{
				has2 = true;
			}

		}

		if (has1 && has2) // eledges[j] is on face
		{
			// fedges.Append (eledges[j]);
			for (int k = 0;k < nfa_ref_edges;k++)
			{
			int w1 = el[ref_faces[fa][fa_ref_edges[k][0] - 1] - 1];
			int w2 = el[ref_faces[fa][fa_ref_edges[k][1] - 1] - 1];

			if (withorientation)
			{
				if (w1 == vi1 && w2 == vi2)
				{
				  fedges[k] = eledges[j];
				}
				if (w1 == vi2 && w2 == vi1)
				{
				  fedges[k] = -eledges[j];
				}
			}
			else
			{
			  if ((w1 == vi1 && w2 == vi2) || (w1 == vi2 && w2 == vi1))
			  {
				fedges[k] = eledges[j];
			  }
			}
			}
		}
		}

		// *testout << " Face " << fnr << endl;
		// *testout << " GetFaceEdges " << fedges << endl;

		return;
	}
	}

	int surfel = GetFace2SurfaceElement(fnr);
	if (surfel != 0)
	{
	GetSurfaceElementEdges(surfel, fedges);
	return;
	}
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: ELEMENT_TYPE GetFaceType(int fnr) const
  public ELEMENT_TYPE GetFaceType(int fnr)
  {
	  return (face2vert.Get fnr[3] == 0) ? ELEMENT_TYPE.TRIG : ELEMENT_TYPE.QUAD;
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetSurfaceElementEdges(int elnr, Array<int> & eledges) const
  public void GetSurfaceElementEdges(int elnr, Array<int> eledges)
  {
	int ned = GetNEdges(mesh.SurfaceElement(elnr).GetType());
	eledges.SetSize(ned);
	for (int i = 0; i < ned; i++)
	{
	  // eledges[i] = abs (surfedges.Get(elnr)[i]);
	  eledges[i] = surfedges.Get(elnr)[i].nr + 1;
	}
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetSurfaceElementFace(int elnr) const
  public int GetSurfaceElementFace(int elnr)
  {
	return surffaces.Get(elnr).fnr + 1;
  }


  /*
  int MeshTopology :: GetFace (SurfaceElementIndex elnr) const
  {
    return surffaces[elnr].fnr;
  }
  */


//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetSurfaceElementEdgeOrientations(int elnr, Array<int> & eorient) const
  public void GetSurfaceElementEdgeOrientations(int elnr, Array<int> eorient)
  {
	int ned = GetNEdges(mesh.SurfaceElement(elnr).GetType());
	eorient.SetSize(ned);
	for (int i = 0; i < ned; i++)
	{
	  // eorient[i] = (surfedges.Get(elnr)[i] > 0) ? 1 : -1;
	  // eorient[i] = (surfedges.Get(elnr)[i].orient) ? -1 : 1;
	  eorient[i] = GetSurfaceElementEdgeOrientation(elnr, i) != 0 ? -1 : 1;
	}
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetSurfaceElementFaceOrientation(int elnr) const
  public int GetSurfaceElementFaceOrientation(int elnr)
  {
	// return (surffaces.Get(elnr)-1) % 8;
	// return surffaces.Get(elnr).forient;
	return GetSurfaceElementFaceOrientation2(elnr);
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetEdges(SurfaceElementIndex elnr, Array<int> & eledges) const
  public void GetEdges(SurfaceElementIndex elnr, Array<int> eledges)
  {
	int ned = GetNEdges(mesh[elnr].GetType());
	eledges.SetSize(ned);
	for (int i = 0; i < ned; i++)
	{
	  // eledges[i] = abs (surfedges[elnr][i])-1;
	  eledges[i] = surfedges[elnr][i].nr;
	}
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetFace(SurfaceElementIndex elnr) const
  public int GetFace(SurfaceElementIndex elnr)
  {
	  return surffaces[elnr].fnr;
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetSurfaceElementEdges(int elnr, int * eledges, int * orient) const
  public int GetSurfaceElementEdges(int elnr, int[] eledges, int[] orient)
  {
	int i;
	if (mesh.GetDimension() == 3 || true)
	{
	if (orient != 0)
	{
		for (i = 0; i < 4; i++)
		{
				/*
		if (!surfedges.Get(elnr)[i]) return i;
		eledges[i] = abs (surfedges.Get(elnr)[i]);
		orient[i] = (surfedges.Get(elnr)[i] > 0 ) ? 1 : -1;
				*/
		if (surfedges.Get(elnr)[i].nr == -1)
		{
			return i;
		}
		eledges[i] = surfedges.Get(elnr)[i].nr + 1;
		// orient[i] = (surfedges.Get(elnr)[i].orient) ? -1 : 1;
				orient[i] = GetSurfaceElementEdgeOrientation(elnr, i) != 0 ? -1 : 1;

		}
	}
	else
	{
		for (i = 0; i < 4; i++)
		{
				/*
		if (!surfedges.Get(elnr)[i]) return i;
		eledges[i] = abs (surfedges.Get(elnr)[i]);
				*/
		if (surfedges.Get(elnr)[i].nr == -1)
		{
			return i;
		}
		eledges[i] = surfedges.Get(elnr)[i].nr + 1;
		}
	}
	return 4;
	}
	else
	{
		/*
	eledges[0] = abs (segedges.Get(elnr));
	if (orient)
	  orient[0] = segedges.Get(elnr) > 0 ? 1 : -1;
		*/
	eledges[0] = segedges.Get(elnr).nr + 1;
	if (orient != 0)
	{
	  // orient[0] = segedges.Get(elnr).orient ? -1 : 1;
		  orient[0] = GetSegmentEdgeOrientation(elnr) != 0 ? -1 : 1;
	}
	}
	return 1;
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const T_EDGE * GetElementEdgesPtr(int elnr) const
  public T_EDGE GetElementEdgesPtr(int elnr)
  {
	  return edges[elnr][0];
  }
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const T_EDGE * GetSurfaceElementEdgesPtr(int selnr) const
  public T_EDGE GetSurfaceElementEdgesPtr(int selnr)
  {
	  return surfedges[selnr][0];
  }
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const T_EDGE * GetSegmentElementEdgesPtr(int selnr) const
  public T_EDGE GetSegmentElementEdgesPtr(int selnr)
  {
	  return segedges[selnr];
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const T_FACE * GetElementFacesPtr(int elnr) const
  public T_FACE GetElementFacesPtr(int elnr)
  {
	  return faces[elnr][0];
  }
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: const T_FACE * GetSurfaceElementFacesPtr(int selnr) const
  public T_FACE GetSurfaceElementFacesPtr(int selnr)
  {
	  return surffaces[selnr];
  }


//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetSurface2VolumeElement(int selnr, int & elnr1, int & elnr2) const
  public void GetSurface2VolumeElement(int selnr, ref int elnr1, ref int elnr2)
  {
	elnr1 = surf2volelement.Get(selnr)[0];
	elnr2 = surf2volelement.Get(selnr)[1];
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetFace2SurfaceElement(int fnr) const
  public int GetFace2SurfaceElement(int fnr)
  {
	  return face2surfel[fnr - 1];
  }


  /*
  ELEMENT_TYPE MeshTopology :: GetFaceType (int fnr) const
  {
    if (face2vert.Get(fnr)[3] == 0) return TRIG; else return QUAD;
  }
  */


//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetVertexElements(int vnr, Array<ElementIndex> & elements) const
  public void GetVertexElements(int vnr, Array<ElementIndex> elements)
  {
	if (vert2element.Size() != 0)
	{
	int ne = vert2element.EntrySize(vnr);
	elements.SetSize(ne);
	for (int i = 1; i <= ne; i++)
	{
	  elements.Elem(i) = vert2element.Get(vnr, i);
	}
	}
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: FlatArray<ElementIndex> GetVertexElements(int vnr) const
  public FlatArray<ElementIndex> GetVertexElements(int vnr)
  {
	  return vert2element[vnr];
  }


  /*
  FlatArray<ElementIndex> MeshTopology :: GetVertexElements (int vnr) const
  {
    if (vert2element)
      return (*vert2element)[vnr];
    return FlatArray<ElementIndex> (0,0);
  }

  FlatArray<SurfaceElementIndex> MeshTopology :: GetVertexSurfaceElements (int vnr) const
  {
    if (vert2surfelement)
      return (*vert2surfelement)[vnr];
    return FlatArray<SurfaceElementIndex> (0,0);
  }

  FlatArray<SegmentIndex> MeshTopology :: GetVertexSegments (int vnr) const
  {
    if (vert2segment)
      return (*vert2segment)[vnr];
    return FlatArray<SegmentIndex> (0,0);
  }
  */

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetVertexSurfaceElements(int vnr, Array<SurfaceElementIndex> & elements) const
  public void GetVertexSurfaceElements(int vnr, Array<SurfaceElementIndex> elements)
  {
	if (vert2surfelement.Size() != 0)
	{
	int i;
	int ne = vert2surfelement.EntrySize(vnr);
	elements.SetSize(ne);
	for (i = 1; i <= ne; i++)
	{
	  elements.Elem(i) = vert2surfelement.Get(vnr, i);
	}
	}
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: FlatArray<SurfaceElementIndex> GetVertexSurfaceElements(int vnr) const
  public FlatArray<SurfaceElementIndex> GetVertexSurfaceElements(int vnr)
  {
	  return vert2surfelement[vnr];
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: FlatArray<SegmentIndex> GetVertexSegments(int vnr) const
  public FlatArray<SegmentIndex> GetVertexSegments(int vnr)
  {
	  return vert2segment[vnr];
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: FlatArray<int> GetVertexPointElements(int vnr) const
  public FlatArray<int> GetVertexPointElements(int vnr)
  {
	  return vert2pointelement[vnr];
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetVerticesEdge(int v1, int v2) const
  public int GetVerticesEdge(int v1, int v2)
  {
	Array<ElementIndex> elements_v1 = new Array<ElementIndex>();
	Array<int> elementedges = new Array<int>();
	GetVertexElements(v1, elements_v1);
	int edv1;
	int edv2;

	for (int i = 0; i < elements_v1.Size(); i++)
	{
	GetElementEdges(elements_v1[i] + 1, elementedges);
	for (int ed = 0; ed < elementedges.Size(); ed++)
	{
		GetEdgeVertices(elementedges[ed], edv1, edv2);
		if ((edv1 == v1 && edv2 == v2) || (edv1 == v2 && edv2 == v1))
		{
		  return elementedges[ed];
		}
	}
	}

	return -1;
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetSegmentVolumeElements(int segnr, Array<ElementIndex> & volels) const
  public void GetSegmentVolumeElements(int segnr, Array<ElementIndex> volels)
  {
	int v1;
	int v2;
	GetEdgeVertices(GetSegmentEdge(segnr), ref v1, ref v2);
	Array<ElementIndex> volels1 = new Array<ElementIndex>();
	Array<ElementIndex> volels2 = new Array<ElementIndex>();
	GetVertexElements(v1, volels1);
	GetVertexElements(v2, volels2);
	volels.SetSize(0);

	for (int eli1 = 1; eli1 <= volels1.Size(); eli1++)
	{
	  if (volels2.Contains(volels1.Elem(eli1)))
	  {
	volels.Append(volels1.Elem(eli1));
	  }
	}
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetSegmentSurfaceElements(int segnr, Array<SurfaceElementIndex> & els) const
  public void GetSegmentSurfaceElements(int segnr, Array<SurfaceElementIndex> els)
  {
	int v1;
	int v2;
	GetEdgeVertices(GetSegmentEdge(segnr), ref v1, ref v2);
	Array<SurfaceElementIndex> els1 = new Array<SurfaceElementIndex>();
	Array<SurfaceElementIndex> els2 = new Array<SurfaceElementIndex>();
	GetVertexSurfaceElements(v1, els1);
	GetVertexSurfaceElements(v2, els2);
	els.SetSize(0);

	for (int eli1 = 1; eli1 <= els1.Size(); eli1++)
	{
	  if (els2.Contains(els1.Elem(eli1)))
	  {
	els.Append(els1.Elem(eli1));
	  }
	}
  }
}

}



