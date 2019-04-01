#if PARALLEL


#define M_PI
#define WIN32_LEAN_AND_MEAN
#define PACKAGE_VERSION
//C++ TO C# CONVERTER TODO TASK: #define macros defined in multiple preprocessor conditionals can only be replaced within the scope of the preprocessor conditional:
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define DLL_HEADER __declspec(dllexport)
	  #define DLL_HEADER
//C++ TO C# CONVERTER TODO TASK: #define macros defined in multiple preprocessor conditionals can only be replaced within the scope of the preprocessor conditional:
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define DLL_HEADER __declspec(dllimport)
	  #define DLL_HEADER
	  #define DLL_HEADER
	  #define DLL_HEADER
//C++ TO C# CONVERTER TODO TASK: #define macros defined in multiple preprocessor conditionals can only be replaced within the scope of the preprocessor conditional:
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define __assume(cond) if (!(cond)) __builtin_unreachable(); else;
#define __assume
//C++ TO C# CONVERTER TODO TASK: #define macros defined in multiple preprocessor conditionals can only be replaced within the scope of the preprocessor conditional:
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define __assume(cond)
#define __assume
//C++ TO C# CONVERTER TODO TASK: #define macros defined in multiple preprocessor conditionals can only be replaced within the scope of the preprocessor conditional:
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define NG_INLINE __forceinline inline
#define NG_INLINE
//C++ TO C# CONVERTER TODO TASK: #define macros defined in multiple preprocessor conditionals can only be replaced within the scope of the preprocessor conditional:
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define NG_INLINE __forceinline inline
#define NG_INLINE
//C++ TO C# CONVERTER TODO TASK: #define macros defined in multiple preprocessor conditionals can only be replaced within the scope of the preprocessor conditional:
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define NG_INLINE inline
#define NG_INLINE
#define VLA
//C++ TO C# CONVERTER TODO TASK: #define macros defined in multiple preprocessor conditionals can only be replaced within the scope of the preprocessor conditional:
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define NG_INLINE inline
#define NG_INLINE
#define noDEMOVERSION
#define noDEVELOP
#define noSTEP
#define noSOLIDGEOM
#define noDEMOAPP
#define noMODELLER
#define noSTAT_STREAM
#define noLOG_STREAM
//C++ TO C# CONVERTER TODO TASK: #define macros defined in multiple preprocessor conditionals can only be replaced within the scope of the preprocessor conditional:
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define NGCORE_API_EXPORT __declspec(dllexport)
		#define NGCORE_API_EXPORT
//C++ TO C# CONVERTER TODO TASK: #define macros defined in multiple preprocessor conditionals can only be replaced within the scope of the preprocessor conditional:
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define NGCORE_API_IMPORT __declspec(dllimport)
		#define NGCORE_API_IMPORT
		#define NGCORE_API_EXPORT
		#define NGCORE_API_IMPORT
//C++ TO C# CONVERTER TODO TASK: #define macros defined in multiple preprocessor conditionals can only be replaced within the scope of the preprocessor conditional:
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define NGCORE_API NGCORE_API_EXPORT
		#define NGCORE_API
//C++ TO C# CONVERTER TODO TASK: #define macros defined in multiple preprocessor conditionals can only be replaced within the scope of the preprocessor conditional:
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define NGCORE_API NGCORE_API_IMPORT
		#define NGCORE_API
//C++ TO C# CONVERTER TODO TASK: #define macros defined in multiple preprocessor conditionals can only be replaced within the scope of the preprocessor conditional:
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define NETGEN_INLINE __forceinline inline
	#define NETGEN_INLINE
	#define NETGEN_LAMBDA_INLINE
//C++ TO C# CONVERTER TODO TASK: #define macros defined in multiple preprocessor conditionals can only be replaced within the scope of the preprocessor conditional:
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define NETGEN_INLINE __forceinline inline
	#define NETGEN_INLINE
	#define NETGEN_LAMBDA_INLINE
//C++ TO C# CONVERTER TODO TASK: #define macros defined in multiple preprocessor conditionals can only be replaced within the scope of the preprocessor conditional:
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define NETGEN_INLINE inline
	#define NETGEN_INLINE
	#define NETGEN_LAMBDA_INLINE
	#define NETGEN_VLA
//C++ TO C# CONVERTER TODO TASK: #define macros defined in multiple preprocessor conditionals can only be replaced within the scope of the preprocessor conditional:
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define NETGEN_INLINE inline
	#define NETGEN_INLINE
	#define NETGEN_LAMBDA_INLINE
//C++ TO C# CONVERTER TODO TASK: #define macros defined in multiple preprocessor conditionals can only be replaced within the scope of the preprocessor conditional:
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define NETGEN_CORE_NGEXEPTION_STR_HELPER(x) #x
#define NETGEN_CORE_NGEXEPTION_STR_HELPER
//C++ TO C# CONVERTER TODO TASK: #define macros defined in multiple preprocessor conditionals can only be replaced within the scope of the preprocessor conditional:
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define NETGEN_CORE_NGEXEPTION_STR(x) NETGEN_CORE_NGEXEPTION_STR_HELPER(x)
#define NETGEN_CORE_NGEXEPTION_STR
//C++ TO C# CONVERTER TODO TASK: #define macros defined in multiple preprocessor conditionals can only be replaced within the scope of the preprocessor conditional:
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define NG_EXCEPTION(s) ngcore::Exception(__FILE__ ":" NETGEN_CORE_NGEXEPTION_STR(__LINE__) "\t"+std::string(s))
#define NG_EXCEPTION
//C++ TO C# CONVERTER TODO TASK: #define macros defined in multiple preprocessor conditionals can only be replaced within the scope of the preprocessor conditional:
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define NETGEN_CHECK_RANGE(value, min, max) { if ((value)<(min) || (value)>=(max)) throw ngcore::RangeException(__FILE__ ":" NETGEN_CORE_NGEXEPTION_STR(__LINE__) "\t", (value), (min), (max)); }
#define NETGEN_CHECK_RANGE
//C++ TO C# CONVERTER TODO TASK: #define macros defined in multiple preprocessor conditionals can only be replaced within the scope of the preprocessor conditional:
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define NETGEN_CHECK_RANGE(value, min, max)
#define NETGEN_CHECK_RANGE
#define SPDLOG_DEBUG_ON
//C++ TO C# CONVERTER TODO TASK: #define macros defined in multiple preprocessor conditionals can only be replaced within the scope of the preprocessor conditional:
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define NETGEN_DEBUG_LOG(logger, ...) SPDLOG_DEBUG(logger, __VA_ARGS__)
#define NETGEN_DEBUG_LOG
//C++ TO C# CONVERTER TODO TASK: #define macros defined in multiple preprocessor conditionals can only be replaced within the scope of the preprocessor conditional:
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define NETGEN_DEBUG_LOG(logger, ...)
#define NETGEN_DEBUG_LOG
#define OMPI_SKIP_MPICXX
#define CLOCKS_PER_SEC
//C++ TO C# CONVERTER TODO TASK: #define macros defined in multiple preprocessor conditionals can only be replaced within the scope of the preprocessor conditional:
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define VT_USER_START(n)
  #define VT_USER_START
//C++ TO C# CONVERTER TODO TASK: #define macros defined in multiple preprocessor conditionals can only be replaced within the scope of the preprocessor conditional:
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define VT_USER_END(n)
  #define VT_USER_END
//C++ TO C# CONVERTER TODO TASK: #define macros defined in multiple preprocessor conditionals can only be replaced within the scope of the preprocessor conditional:
//C++ TO C# CONVERTER NOTE: The following #define macro was replaced in-line:
//ORIGINAL LINE: #define VT_TRACER(n)
  #define VT_TRACER
#define GZSTREAM_H
#define EPSGEOM
#define ELEMENT_MAXPOINTS
#define ELEMENT2D_MAXPOINTS
#define MULTIPOINTGEOMINFO_MAX
#define _INCLUDE_MORE

namespace netgen
{


  public class ParallelMeshTopology : System.IDisposable
  {
	private readonly Mesh mesh;

	/**
	   mapping from local to distant vertex number
	   each row of the table corresponds to one vertex
	   each row contains a list of pairs (procnr, dist_vnum)
	*/

	private TABLE<int> loc2distvert = new TABLE<int>();
	private TABLE<int> loc2distedge = new TABLE<int>();
	private TABLE<int> loc2distface = new TABLE<int>();

	private Array<int> glob_vert = new Array<int>();
	private Array<int> glob_edge = new Array<int>();
	private Array<int> glob_face = new Array<int>();
	private Array<int> glob_el = new Array<int>();
	private Array<int> glob_surfel = new Array<int>();
	private Array<int> glob_segm = new Array<int>();

	private bool is_updated;


	public ParallelMeshTopology(Mesh amesh)
	{
		this.mesh = new netgen.Mesh(amesh);
	  is_updated = false;
	}

	public void Dispose()
	{
	  ;
	}

	public void Reset()
	{
	  *testout << "ParallelMeshTopology::Reset" << "\n";

	  NgMPI_Comm comm = mesh.GetCommunicator();
	  int id = comm.Rank();
	  int ntasks = comm.Size();

	  if (ntasks == 1)
	  {
		  return;
	  }

	  int ned = mesh.GetTopology().GetNEdges();
	  int nfa = mesh.GetTopology().GetNFaces();

	  if (glob_edge.Size() != ned)
	  {
	  glob_edge.SetSize(ned);
	  glob_face.SetSize(nfa);
	  glob_edge = -1;
	  glob_face = -1;

	  loc2distedge.ChangeSize(ned);
	  loc2distface.ChangeSize(nfa);
	  }

	  if (glob_vert.Size() != mesh.GetNV())
	  {
	  SetNV(mesh.GetNV());
	  SetNE(mesh.GetNE());
	  }
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void Print() const
	public void Print()
	{
	  ;
	}

//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	private int UpdateCoarseGrid_timer = NgProfiler.CreateTimer("UpdateCoarseGrid");
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	private int UpdateCoarseGrid_timere = NgProfiler.CreateTimer("UpdateCoarseGrid - ex edges");
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
	private int UpdateCoarseGrid_timerf = NgProfiler.CreateTimer("UpdateCoarseGrid - ex faces");

	public void UpdateCoarseGrid()
	{
	  // cout << "UpdateCoarseGrid" << endl;
	  // if (is_updated) return;

	  NgMPI_Comm comm = mesh.GetCommunicator();
	  int id = comm.Rank();
	  int ntasks = comm.Size();

	  if (ntasks == 1)
	  {
		  return;
	  }

	  Reset();
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	  static int timer = NgProfiler::CreateTimer("UpdateCoarseGrid");
	  NgProfiler.RegionTimer reg = new NgProfiler.RegionTimer(UpdateCoarseGrid_timer);


	  (*testout) << "UPDATE COARSE GRID PARALLEL TOPOLOGY " << "\n";
	  if (id == 0)
	  {
		PrintMessage(1, "update parallel topology");
	  }


	  // UpdateCoarseGridGlobal();



	  // MPI_Barrier (MPI_COMM_WORLD);

	  MPI_Group MPI_GROUP_comm = new MPI_Group();
	  MPI_Group MPI_LocalGroup = new MPI_Group();
	  MPI_Comm MPI_LocalComm = new MPI_Comm();

	  int[] process_ranks = {0};
	  MPI_Comm_group(comm, MPI_GROUP_comm);
	  MPI_Group_excl(MPI_GROUP_comm, 1, process_ranks, MPI_LocalGroup);
	  MPI_Comm_create(comm, MPI_LocalGroup, MPI_LocalComm);

	  if (id == 0)
	  {
		  return;
	  }

	  MeshTopology topology = mesh.GetTopology();

	  Array<int> cnt_send = new Array<int>(ntasks - 1);


	  // update new vertices after mesh-refinement
	  if (mesh.mlbetweennodes.Size() > 0)
	  {
	  // cout << "UpdateCoarseGrid - vertices" << endl;
		  int newnv = mesh.mlbetweennodes.Size();
		  loc2distvert.ChangeSize(mesh.mlbetweennodes.Size());
	  /*
		  for (PointIndex pi = PointIndex::BASE; pi < newnv+PointIndex::BASE; pi++)
			{
			  PointIndex v1 = mesh.mlbetweennodes[pi][0];
			  PointIndex v2 = mesh.mlbetweennodes[pi][1];
			  if (mesh.mlbetweennodes[pi][0] != PointIndex::BASE-1)
				for (int dest = 1; dest < ntasks; dest++)
				  if (IsExchangeVert (dest, v1) && IsExchangeVert (dest, v2))
					SetDistantPNum(dest, pi);
			}
	  */

	  bool changed = true;
	  while (changed)
	  {
		  changed = false;

		  // build exchange vertices
		  cnt_send = 0;
		  foreach (PointIndex pi in mesh.Points().Range())
		  {
			foreach (int dist in GetDistantPNums(pi - PointIndex.BASE))
			{
		  cnt_send[dist - 1]++;
			}
		  }
		  TABLE<int> dest2vert = new TABLE<int>(cnt_send);
		  foreach (PointIndex pi in mesh.Points().Range())
		  {
			foreach (int dist in GetDistantPNums(pi - PointIndex.BASE))
			{
		  dest2vert.Add(dist - 1, pi);
			}
		  }


		  for (PointIndex pi = PointIndex.BASE; pi < newnv + PointIndex.BASE; pi++)
		  {
		  PointIndex v1 = mesh.mlbetweennodes[pi][0];
		  PointIndex v2 = mesh.mlbetweennodes[pi][1];
		  if (mesh.mlbetweennodes[pi][0] != PointIndex.BASE-1)
		  {
			// for (int dest = 1; dest < ntasks; dest++)
			foreach (int dest in GetDistantPNums(v1 - PointIndex.BASE))
			{
			  if (IsExchangeVert(dest, new netgen.PointIndex(v1)) && IsExchangeVert(dest, new netgen.PointIndex(v2)))
			  {
				cnt_send[dest - 1]++;
			  }
			}
		  }
		  }

		  TABLE<int> dest2pair = new TABLE<int>(cnt_send);
		  // for (int dest = 1; dest < ntasks; dest++)
			for (PointIndex pi = PointIndex.BASE; pi < newnv + PointIndex.BASE; pi++)
			{
			PointIndex v1 = mesh.mlbetweennodes[pi][0];
			PointIndex v2 = mesh.mlbetweennodes[pi][1];
			if (mesh.mlbetweennodes[pi][0] != PointIndex.BASE-1)
			{
			  foreach (int dest in GetDistantPNums(v1 - PointIndex.BASE))
			  {
				if (IsExchangeVert(dest, new netgen.PointIndex(v1)) && IsExchangeVert(dest, new netgen.PointIndex(v2)))
				{
			  dest2pair.Add(dest - 1, pi);
				}
			  }
			}
			}

		  cnt_send = 0;
		  int v1;
		  int v2;
		  for (PointIndex pi = PointIndex.BASE; pi < newnv + PointIndex.BASE; pi++)
		  {
		  PointIndex v1 = mesh.mlbetweennodes[pi][0];
		  PointIndex v2 = mesh.mlbetweennodes[pi][1];
		  if (mesh.mlbetweennodes[pi][0] != PointIndex.BASE-1)
		  {
			foreach (int dest in GetDistantPNums(v1 - PointIndex.BASE))
			{
			  if (IsExchangeVert(dest, new netgen.PointIndex(v2)))
			  {
				cnt_send[dest - 1] += 2;
			  }
			}
		  }
		  }

		  TABLE<int> send_verts = new TABLE<int>(cnt_send);

		  Array<int, PointIndex.BASE> loc2exchange = new Array<int, PointIndex.BASE>(mesh.GetNV());
		  for (int dest = 1; dest < ntasks; dest++)
		  {
			if (dest != id)
			{
			loc2exchange = -1;
			int cnt = 0;
			/*
			for (PointIndex pi : mesh.Points().Range())
			  if (IsExchangeVert(dest, pi))
			    loc2exchange[pi] = cnt++;
			*/
			foreach (PointIndex pi in dest2vert[dest - 1])
			{
			  loc2exchange[pi] = cnt++;
			}

			// for (PointIndex pi = PointIndex::BASE; pi < newnv+PointIndex::BASE; pi++)
			foreach (PointIndex pi in dest2pair[dest - 1])
			{
				PointIndex v1 = mesh.mlbetweennodes[pi][0];
				PointIndex v2 = mesh.mlbetweennodes[pi][1];
				if (mesh.mlbetweennodes[pi][0] != PointIndex.BASE-1)
				{
			  if (IsExchangeVert(dest, new netgen.PointIndex(v1)) && IsExchangeVert(dest, new netgen.PointIndex(v2)))
			  {
				  send_verts.Add(dest - 1, loc2exchange[v1]);
				  send_verts.Add(dest - 1, loc2exchange[v2]);
			  }
				}
			}
			}
		  }

		  TABLE<int> recv_verts = new TABLE<int>(ntasks - 1);
		  netgen.GlobalMembers.MyMPI_ExchangeTable(send_verts, recv_verts, MPI_TAG_MESH + 9, MPI_LocalComm);

		  for (int dest = 1; dest < ntasks; dest++)
		  {
			if (dest != id)
			{
			loc2exchange = -1;
			int cnt = 0;
			/*
			for (PointIndex pi : mesh.Points().Range())
			  if (IsExchangeVert(dest, pi))
			    loc2exchange[pi] = cnt++;
			*/
			foreach (PointIndex pi in dest2vert[dest - 1])
			{
			  loc2exchange[pi] = cnt++;
			}

			FlatArray<int> recvarray = recv_verts[dest - 1];
			for (int ii = 0; ii < recvarray.Size(); ii += 2)
			{
			  foreach (PointIndex pi in dest2pair[dest - 1])
			  {
				// for (PointIndex pi = PointIndex::BASE; pi < newnv+PointIndex::BASE; pi++)
			  PointIndex v1 = mesh.mlbetweennodes[pi][0];
			  PointIndex v2 = mesh.mlbetweennodes[pi][1];
			  if (mesh.mlbetweennodes[pi][0] != PointIndex.BASE-1)
			  {
				  INDEX_2 re = new INDEX_2(recvarray[ii], recvarray[ii + 1]);
				  INDEX_2 es = new INDEX_2(loc2exchange[v1], loc2exchange[v2]);
				  if (es == re && !IsExchangeVert(dest, new netgen.PointIndex(pi)))
				  {
				  SetDistantPNum(dest, new netgen.PointIndex(pi));
				  changed = true;
				  }
			  }
			  }
			}
			}
		  }
	  }
	  }

	  Array<int> sendarray = new Array<int>();
	  Array<int> recvarray = new Array<int>();
	  // cout << "UpdateCoarseGrid - edges" << endl;

	  // static int timerv = NgProfiler::CreateTimer ("UpdateCoarseGrid - ex vertices");
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	  static int timere = NgProfiler::CreateTimer("UpdateCoarseGrid - ex edges");
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	  static int timerf = NgProfiler::CreateTimer("UpdateCoarseGrid - ex faces");


	  NgProfiler.StartTimer(UpdateCoarseGrid_timere);


	  int nfa = topology.GetNFaces();
	  int ned = topology.GetNEdges();

	  // build exchange vertices
	  cnt_send = 0;
	  foreach (PointIndex pi in mesh.Points().Range())
	  {
		foreach (int dist in GetDistantPNums(pi - PointIndex.BASE))
		{
	  cnt_send[dist - 1]++;
		}
	  }
	  TABLE<int> dest2vert = new TABLE<int>(cnt_send);
	  foreach (PointIndex pi in mesh.Points().Range())
	  {
		foreach (int dist in GetDistantPNums(pi - PointIndex.BASE))
		{
	  dest2vert.Add(dist - 1, pi);
		}
	  }

	  // exchange edges
	  cnt_send = 0;
	  int v1;
	  int v2;
	  for (int edge = 1; edge <= ned; edge++)
	  {
	  topology.GetEdgeVertices(edge, ref v1, ref v2);
	  for (int dest = 1; dest < ntasks; dest++)
	  {
		if (IsExchangeVert(dest, v1) && IsExchangeVert(dest, v2))
		{
		  cnt_send[dest - 1] += 1;
		}
	  }
	  }

	  TABLE<int> dest2edge = new TABLE<int>(cnt_send);
	  foreach (int v in cnt_send)
	  {
		  v *= 2;
	  }
	  TABLE<int> send_edges = new TABLE<int>(cnt_send);

	  for (int edge = 1; edge <= ned; edge++)
	  {
	  topology.GetEdgeVertices(edge, ref v1, ref v2);
	  for (int dest = 1; dest < ntasks; dest++)
	  {
		if (IsExchangeVert(dest, v1) && IsExchangeVert(dest, v2))
		{
		  dest2edge.Add(dest - 1, edge);
		}
	  }
	  }


	  Array<int, PointIndex.BASE> loc2exchange = new Array<int, PointIndex.BASE>(mesh.GetNV());
	  for (int dest = 1; dest < ntasks; dest++)
	  {
		  loc2exchange = -1;
		  int cnt = 0;
		  foreach (PointIndex pi in dest2vert[dest - 1])
		  {
		loc2exchange[pi] = cnt++;
		  }

	  foreach (int edge in dest2edge[dest - 1])
	  {
			  topology.GetEdgeVertices(edge, ref v1, ref v2);
			  if (IsExchangeVert(dest, v1) && IsExchangeVert(dest, v2))
			  {
				  send_edges.Add(dest - 1, loc2exchange[v1]);
				  send_edges.Add(dest - 1, loc2exchange[v2]);
			  }
	  }
	  }

	  // cout << "UpdateCoarseGrid - edges mpi-exchange" << endl;
	  TABLE<int> recv_edges = new TABLE<int>(ntasks - 1);
	  netgen.GlobalMembers.MyMPI_ExchangeTable(send_edges, recv_edges, MPI_TAG_MESH + 9, MPI_LocalComm);
	  // cout << "UpdateCoarseGrid - edges mpi-exchange done" << endl;

	  /*
	  for (int dest = 1; dest < ntasks; dest++)
	    {
	  auto ex2loc = dest2vert[dest-1];
	  FlatArray<int> recvarray = recv_edges[dest-1];
	      for (int ii = 0; ii < recvarray.Size(); ii+=2)
		for (int edge : dest2edge[dest-1])
		  {
			topology.GetEdgeVertices (edge, v1, v2);
			INDEX_2 re(ex2loc[recvarray[ii]],
			   ex2loc[recvarray[ii+1]]);
			INDEX_2 es(v1, v2);
			if (es == re)
		  SetDistantEdgeNum(dest, edge);
		  }
	    }
	  */

	  for (int dest = 1; dest < ntasks; dest++)
	  {
	  var ex2loc = dest2vert[dest - 1];
	  if (ex2loc.Size() == 0)
	  {
		  continue;
	  }

	  INDEX_2_CLOSED_HASHTABLE<int> vert2edge = new INDEX_2_CLOSED_HASHTABLE<int>((uint)(2 * dest2edge[dest - 1].Size() + 10));
	  foreach (int edge in dest2edge[dest - 1])
	  {
		  topology.GetEdgeVertices(edge, ref v1, ref v2);
		  vert2edge.Set(new INDEX_2(v1, v2), edge);
	  }

	  FlatArray<int> recvarray = recv_edges[dest - 1];
		  for (int ii = 0; ii < recvarray.Size(); ii += 2)
		  {
		  INDEX_2 re = new INDEX_2(ex2loc[recvarray[ii]], ex2loc[recvarray[ii + 1]]);
		  if (vert2edge.Used(re))
		  {
			SetDistantEdgeNum(dest, vert2edge.Get(re));
		  }
		  }
	  }



	  NgProfiler.StopTimer(UpdateCoarseGrid_timere);

	  // MPI_Barrier (MPI_LocalComm);

	  // cout << "UpdateCoarseGrid - faces" << endl;
	  if (mesh.GetDimension() == 3)
	  {
	  NgProfiler.StartTimer(UpdateCoarseGrid_timerf);
	  Array<int> verts = new Array<int>();

	  // exchange faces
	  cnt_send = 0;
	  for (int face = 1; face <= nfa; face++)
	  {
		  topology.GetFaceVertices(face, verts);
		  for (int dest = 1; dest < ntasks; dest++)
		  {
			if (dest != id)
			{
		  if (IsExchangeVert(dest, verts[0]) && IsExchangeVert(dest, verts[1]) && IsExchangeVert(dest, verts[2]))
		  {
			cnt_send[dest - 1]++;
		  }
			}
		  }
	  }

	  TABLE<int> dest2face = new TABLE<int>(cnt_send);
	  for (int face = 1; face <= nfa; face++)
	  {
		  topology.GetFaceVertices(face, verts);
		  for (int dest = 1; dest < ntasks; dest++)
		  {
			if (dest != id)
			{
		  if (IsExchangeVert(dest, verts[0]) && IsExchangeVert(dest, verts[1]) && IsExchangeVert(dest, verts[2]))
		  {
			dest2face.Add(dest - 1, face);
		  }
			}
		  }
	  }

	  foreach (int c in cnt_send)
	  {
		  c *= 3;
	  }
	  TABLE<int> send_faces = new TABLE<int>(cnt_send);
	  Array<int, PointIndex.BASE> loc2exchange = new Array<int, PointIndex.BASE>(mesh.GetNV());
	  for (int dest = 1; dest < ntasks; dest++)
	  {
		if (dest != id)
		{
			/*
			loc2exchange = -1;
			int cnt = 0;
			for (PointIndex pi : mesh.Points().Range())
		  if (IsExchangeVert(dest, pi))
			loc2exchange[pi] = cnt++;
			*/
			if (dest2vert[dest - 1].Size() == 0)
			{
				continue;
			}

			loc2exchange = -1;
			int cnt = 0;
			foreach (PointIndex pi in dest2vert[dest - 1])
			{
		  loc2exchange[pi] = cnt++;
			}

			foreach (int face in dest2face[dest - 1])
			{
			topology.GetFaceVertices(face, verts);
			if (IsExchangeVert(dest, verts[0]) && IsExchangeVert(dest, verts[1]) && IsExchangeVert(dest, verts[2]))
			{
				send_faces.Add(dest - 1, loc2exchange[verts[0]]);
				send_faces.Add(dest - 1, loc2exchange[verts[1]]);
				send_faces.Add(dest - 1, loc2exchange[verts[2]]);
			}
			}
		}
	  }

	  // cout << "UpdateCoarseGrid - faces mpi-exchange" << endl;
	  TABLE<int> recv_faces = new TABLE<int>(ntasks - 1);
	  netgen.GlobalMembers.MyMPI_ExchangeTable(send_faces, recv_faces, MPI_TAG_MESH + 9, MPI_LocalComm);
	  // cout << "UpdateCoarseGrid - faces mpi-exchange done" << endl;

	  /*
	  for (int dest = 1; dest < ntasks; dest++)
	    if (dest != id)
	      {
	        loc2exchange = -1;
	        int cnt = 0;
	        for (PointIndex pi : dest2vert[dest-1])
	  	loc2exchange[pi] = cnt++;
	  
	        FlatArray<int> recvarray = recv_faces[dest-1];
	        for (int ii = 0; ii < recvarray.Size(); ii+=3)
	  	for (int face : dest2face[dest-1])
	  	  {
	  	    topology.GetFaceVertices (face, verts);
	  	    INDEX_3 re(recvarray[ii], recvarray[ii+1], recvarray[ii+2]);
	  	    INDEX_3 es(loc2exchange[verts[0]], loc2exchange[verts[1]], loc2exchange[verts[2]]);
	  	    if (es == re)
	  	      SetDistantFaceNum(dest, face);
	  	  }
	      }
	  */


	  for (int dest = 1; dest < ntasks; dest++)
	  {
		  var ex2loc = dest2vert[dest - 1];
		  if (ex2loc.Size() == 0)
		  {
			  continue;
		  }

		  INDEX_3_CLOSED_HASHTABLE<int> vert2face = new INDEX_3_CLOSED_HASHTABLE<int>(2 * dest2face[dest - 1].Size() + 10);
		  foreach (int face in dest2face[dest - 1])
		  {
		  topology.GetFaceVertices(face, verts);
		  vert2face.Set(new INDEX_3(verts[0], verts[1], verts[2]), face);
		  }

		  FlatArray<int> recvarray = recv_faces[dest - 1];
		  for (int ii = 0; ii < recvarray.Size(); ii += 3)
		  {
		  INDEX_3 re = new INDEX_3(ex2loc[recvarray[ii]], ex2loc[recvarray[ii + 1]], ex2loc[recvarray[ii + 2]]);
		  if (vert2face.Used(re))
		  {
			SetDistantFaceNum(dest, vert2face.Get(re));
		  }
		  }
	  }






	  /*
	    Array<int,1> glob2loc;
  
	  int maxface = 0;
	  for (int face = 1; face <= nfa; face++)
	    maxface = max (maxface, GetGlobalFaceNum (face));
	  
	  // glob2loc.SetSize (nfaglob);
	  glob2loc.SetSize (maxface);
	  glob2loc = -1;
	  
	  for (int loc = 1; loc <= nfa; loc++)
	    glob2loc[GetGlobalFaceNum(loc)] = loc;
	  
	  cnt_send = 0;
	  Array<int> verts;
	  for (int face = 1; face <= nfa; face++)
	    {
	      topology.GetFaceVertices (face, verts);
	      for (int dest = 1; dest < ntasks; dest++)
	        if (IsExchangeVert (dest, verts[0]) &&
	  	  IsExchangeVert (dest, verts[1]) &&
	  	  IsExchangeVert (dest, verts[2]))
	  	{
	  	  cnt_send[dest-1]+=2;
	  	}
	    }
	  
	  TABLE<int> send_faces(cnt_send);
	  for (int face = 1; face <= nfa; face++)
	    {
	      topology.GetFaceVertices (face, verts);
	      for (int dest = 1; dest < ntasks; dest++)
	        {
	  	if (IsExchangeVert (dest, verts[0]) &&
	  	    IsExchangeVert (dest, verts[1]) &&
	  	    IsExchangeVert (dest, verts[2]))
	  	  {
	  	    send_faces.Add (dest-1, GetGlobalFaceNum(face));
	  	    send_faces.Add (dest-1, face);
	  	  }
	        }
	    }
	  TABLE<int> recv_faces(ntasks-1);
	  MyMPI_ExchangeTable (send_faces, recv_faces, MPI_TAG_MESH+8, MPI_LocalComm);
	  
	  for (int sender = 1; sender < ntasks; sender ++)
	    if (id != sender)
	      {
	        FlatArray<int> recvarray = recv_faces[sender-1];
	  
	        for (int ii = 0; ii < recvarray.Size(); )
	  	{
	  	  int globf = recvarray[ii++];
	  	  int distf = recvarray[ii++];
	  
	  	  if (globf <= maxface)
	  	    {
	  	      int locf = glob2loc[globf];
	  	      if (locf != -1)
	  		SetDistantFaceNum (sender, locf);
	  	    }
	  	}
	      }
	  */

	  NgProfiler.StopTimer(UpdateCoarseGrid_timerf);
	  }
	  // cout << "UpdateCoarseGrid - done" << endl;

	  is_updated = true;
	}

	public void UpdateCoarseGridGlobal()
	{
	  // cout << "updatecoarsegridglobal called" << endl;
	  if (id == 0)
	  {
		PrintMessage(3, "UPDATE GLOBAL COARSEGRID STARTS");
	  }

	  int timer = NgProfiler.CreateTimer("UpdateCoarseGridGlobal");
	  NgProfiler.RegionTimer reg = new NgProfiler.RegionTimer(timer);

	  *testout << "ParallelMeshTopology :: UpdateCoarseGridGlobal" << "\n";

	  MeshTopology topology = mesh.GetTopology();
	  MPI_Comm comm = mesh.GetCommunicator();

	  if (id == 0)
	  {
	  Array<Array<int>> sendarrays = new Array<Array<int>>(ntasks);
	  for (int dest = 1; dest < ntasks; dest++)
	  {
		sendarrays[dest] = new Array<int>();
	  }

	  Array<int> edges = new Array<int>();
	  Array<int> faces = new Array<int>();
	  for (int el = 1; el <= mesh.GetNE(); el++)
	  {
		  topology.GetElementFaces(el, faces);
		  topology.GetElementEdges(el, edges);
		  Element volel = mesh.VolumeElement(el);

		  // Array<int> & sendarray = *sendarrays[volel.GetPartition()];
			  Array<int> sendarray = *sendarrays[mesh.vol_partition[el - 1]];

		  for (int i = 0; i < edges.Size(); i++)
		  {
			sendarray.Append(edges[i]);
		  }
		  for (int i = 0; i < faces.Size(); i++)
		  {
			sendarray.Append(faces[i]);
		  }
	  }

	  for (int el = 1; el <= mesh.GetNSE(); el++)
	  {
		  topology.GetSurfaceElementEdges(el, edges);
		  Element2d surfel = mesh.SurfaceElement(el);
		  // Array<int> & sendarray = *sendarrays[surfel.GetPartition()];
			  Array<int> sendarray = *sendarrays[mesh.surf_partition[el - 1]];

		  for (int i = 0; i < edges.Size(); i++)
		  {
			sendarray.Append(edges[i]);
		  }
		  sendarray.Append(topology.GetSurfaceElementFace(el));
	  }

	  Array<MPI_Request> sendrequests = new Array<MPI_Request>();
	  for (int dest = 1; dest < ntasks; dest++)
	  {
		sendrequests.Append(netgen.GlobalMembers.MyMPI_ISend(sendarrays[dest], dest, MPI_TAG_MESH + 10, new MPI_Comm(comm)));
	  }
	  MPI_Waitall(sendrequests.Size(), sendrequests[0], MPI_STATUS_IGNORE);

	  for (int dest = 1; dest < ntasks; dest++)
	  {
		sendarrays[dest] = null;
	  }
	  }

	  else

	  {
	  Array<int> recvarray = new Array<int>();
	  netgen.GlobalMembers.MyMPI_Recv(ref recvarray, 0, MPI_TAG_MESH + 10, new MPI_Comm(comm));

	  int ii = 0;

	  Array<int> faces = new Array<int>();
	  Array<int> edges = new Array<int>();

	  for (int volel = 1; volel <= mesh.GetNE(); volel++)
	  {
		  topology.GetElementEdges(volel, edges);
		  for (int i = 0; i < edges.Size(); i++)
		  {
			SetLoc2Glob_Edge(edges[i], recvarray[ii++]);
		  }

		  topology.GetElementFaces(volel, faces);
		  for (int i = 0; i < faces.Size(); i++)
		  {
			SetLoc2Glob_Face(faces[i], recvarray[ii++]);
		  }
	  }

	  for (int surfel = 1; surfel <= mesh.GetNSE(); surfel++)
	  {
		  topology.GetSurfaceElementEdges(surfel, edges);
		  for (int i = 0; i < edges.Size(); i++)
		  {
			SetLoc2Glob_Edge(edges[i], recvarray[ii++]);
		  }
		  int face = topology.GetSurfaceElementFace(surfel);
		  SetLoc2Glob_Face(face, recvarray[ii++]);
	  }
	  }

	  is_updated = true;
	}

	// bool DoCoarseUpdate() const { return !coarseupdate; }



	/// set number of local vertices, reset sizes of loc2dist_vert, isexchangevert...
	public void SetNV(int anv)
	{
	  glob_vert.SetSize(anv);
	  glob_vert = -1;
	  loc2distvert.ChangeSize(anv);
	}

	public void SetNE(int ane)
	{
	  glob_el.SetSize(ane);
	  glob_el = -1;
	}

	public void SetNSE(int anse)
	{
	  glob_surfel.SetSize(anse);
	  glob_surfel = -1;
	}

	public void SetNSegm(int anseg)
	{
	  glob_segm.SetSize(anseg);
	  glob_segm = -1;
	}


	public void SetLoc2Glob_Vert(int locnum, int globnum)
	{
		glob_vert[locnum - 1] = globnum;
	}
	public void SetLoc2Glob_Edge(int locnum, int globnum)
	{
		glob_edge[locnum - 1] = globnum;
	}
	public void SetLoc2Glob_Face(int locnum, int globnum)
	{
		glob_face[locnum - 1] = globnum;
	}
	public void SetLoc2Glob_VolEl(int locnum, int globnum)
	{
		glob_el[locnum - 1] = globnum;
	}
	public void SetLoc2Glob_SurfEl(int locnum, int globnum)
	{
		glob_surfel[locnum - 1] = globnum;
	}
	public void SetLoc2Glob_Segm(int locnum, int globnum)
	{
		glob_segm[locnum - 1] = globnum;
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetGlobalPNum(int locnum) const
	public int GetGlobalPNum(int locnum)
	{
		return glob_vert[locnum - 1];
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetGlobalEdgeNum(int locnum) const
	public int GetGlobalEdgeNum(int locnum)
	{
		return glob_edge[locnum - 1];
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetGlobalFaceNum(int locnum) const
	public int GetGlobalFaceNum(int locnum)
	{
		return glob_face[locnum - 1];
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetGlobalElNum(int locnum) const
	public int GetGlobalElNum(int locnum)
	{
		return glob_el[locnum - 1];
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetGlobalSElNum(int locnum) const
	public int GetGlobalSElNum(int locnum)
	{
		return glob_surfel[locnum - 1];
	}


	public void SetDistantFaceNum(int dest, int locnum)
	{
	  for (int i = 0; i < loc2distface[locnum - 1].Size(); i += 1)
	  {
		if (loc2distface[locnum - 1][i] == dest)
		{
	  return;
		}
	  }
	  loc2distface.Add(locnum - 1, dest);
	}

	public void SetDistantPNum(int dest, int locnum)
	{
	  for (int i = 0; i < loc2distvert[locnum - 1].Size(); i += 1)
	  {
		if (loc2distvert[locnum - 1][i] == dest)
		{
	  return;
		}
	  }
	  loc2distvert.Add(locnum - 1, dest);
	}

	public void SetDistantEdgeNum(int dest, int locnum)
	{
	  for (int i = 0; i < loc2distedge[locnum - 1].Size(); i += 1)
	  {
		if (loc2distedge[locnum - 1][i] == dest)
		{
	  return;
		}
	  }
	  loc2distedge.Add(locnum - 1, dest);
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetNDistantPNums(int locpnum) const
	public int GetNDistantPNums(int locpnum)
	{
		return loc2distvert[locpnum - 1].Size();
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetNDistantFaceNums(int locfacenum) const
	public int GetNDistantFaceNums(int locfacenum)
	{
		return loc2distface[locfacenum - 1].Size();
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetNDistantEdgeNums(int locedgenum) const
	public int GetNDistantEdgeNums(int locedgenum)
	{
		return loc2distedge[locedgenum - 1].Size();
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetDistantPNums(int locpnum, int * distpnums) const
	public void GetDistantPNums(int locpnum, int[] distpnums)
	{
	  for (int i = 0; i < loc2distvert[locpnum - 1].Size(); i++)
	  {
	distpnums[i] = loc2distvert[locpnum - 1][i];
	  }
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetDistantFaceNums(int locfacenum, int * distfacenums) const
	public void GetDistantFaceNums(int locfacenum, int[] distfacenums)
	{
	  for (int i = 0; i < loc2distface[locfacenum - 1].Size(); i++)
	  {
	distfacenums[i] = loc2distface[locfacenum - 1][i];
	  }
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetDistantFaceNums(int locfacenum, Array<int> & distfacenums) const
	public void GetDistantFaceNums(int locfacenum, ref Array<int> distfacenums)
	{
	  distfacenums = loc2distface[locfacenum - 1];
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetDistantEdgeNums(int locedgenum, int * distedgenums) const
	public void GetDistantEdgeNums(int locedgenum, int[] distedgenums)
	{
	  for (int i = 0; i < loc2distedge[locedgenum - 1].Size(); i++)
	  {
	distedgenums[i] = loc2distedge[locedgenum - 1][i];
	  }
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetDistantEdgeNums(int locedgenum, Array<int> & distedgenums) const
	public void GetDistantEdgeNums(int locedgenum, ref Array<int> distedgenums)
	{
	  distedgenums = loc2distedge[locedgenum - 1];
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: FlatArray<int> GetDistantPNums(int locnum) const
	public FlatArray<int> GetDistantPNums(int locnum)
	{
		return loc2distvert[locnum];
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: FlatArray<int> GetDistantFaceNums(int locnum) const
	public FlatArray<int> GetDistantFaceNums(int locnum)
	{
		return loc2distface[locnum];
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: FlatArray<int> GetDistantEdgeNums(int locnum) const
	public FlatArray<int> GetDistantEdgeNums(int locnum)
	{
		return loc2distedge[locnum];
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: bool IsExchangeVert(int dest, int vnum) const
	public bool IsExchangeVert(int dest, int vnum)
	{
	  return loc2distvert[vnum - 1].Contains(dest);
	}
  }


}









#endif
