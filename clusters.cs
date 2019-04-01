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
/* File:   clusers.hh                                                     */
/* Author: Joachim Schoeberl                                              */
/* Date:   28. Apr. 01                                                    */
/**************************************************************************/

/*
  Anisotropic clusters

  nodes, edges, faces, elements
*/


public class AnisotropicClusters
{
  private readonly Mesh mesh;

  private int nv;
  private int ned;
  private int nfa;
  private int ne;

  // connected nodes, nodes = vertices, edges, faces, elements
  private Array<int> cluster_reps = new Array<int>();

  public AnisotropicClusters(Mesh amesh)
  {
	  this.mesh = amesh;
	;
  }

  public void Dispose()
  {
	;
  }

//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private int Update_timer = NgProfiler.CreateTimer("clusters");
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private int Update_timer1 = NgProfiler.CreateTimer("clusters1");
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private int Update_timer2 = NgProfiler.CreateTimer("clusters2");
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private int Update_timer3 = NgProfiler.CreateTimer("clusters3");

  public void Update(TaskManager tm = DummyTaskManager, Tracer tracer = DummyTracer)
  {
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static int timer = NgProfiler::CreateTimer("clusters");
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static int timer1 = NgProfiler::CreateTimer("clusters1");
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static int timer2 = NgProfiler::CreateTimer("clusters2");
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static int timer3 = NgProfiler::CreateTimer("clusters3");
	NgProfiler.RegionTimer reg = new NgProfiler.RegionTimer(Update_timer);

	MeshTopology top = mesh.GetTopology();

	var id = this.mesh.GetCommunicator().Rank();
	var ntasks = this.mesh.GetCommunicator().Size();

	bool hasedges = top.HasEdges();
	bool hasfaces = top.HasFaces();

	if (!hasedges || !hasfaces)
	{
		return;
	}

	if (id == 0)
	{
	  PrintMessage(3, "Update clusters");
	}

	nv = mesh.GetNV();
	ned = top.GetNEdges();
	nfa = top.GetNFaces();
	ne = mesh.GetNE();
	int nse = mesh.GetNSE();

	cluster_reps.SetSize(nv + ned + nfa + ne);
	cluster_reps = -1;

	Array<int> llist = new Array<int>(nv + ned + nfa + ne);
	llist = 0;

	Array<int> nnums = new Array<int>();
	Array<int> ednums = new Array<int>();
	Array<int> fanums = new Array<int>();
	int changed;

	NgProfiler.StartTimer(Update_timer1);


	/*
	for (int i = 1; i <= ne; i++)
	  {
	const Element & el = mesh.VolumeElement(i);
	ELEMENT_TYPE typ = el.GetType();
	
	top.GetElementEdges (i, ednums);
	top.GetElementFaces (i, fanums);
	
	int elnv = top.GetNVertices (typ);
	int elned = ednums.Size();
	int elnfa = fanums.Size();

	nnums.SetSize(elnv+elned+elnfa+1);
	for (int j = 1; j <= elnv; j++)
	  nnums.Elem(j) = el.PNum(j);
	for (int j = 1; j <= elned; j++)
	  nnums.Elem(elnv+j) = nv+ednums.Elem(j);
	for (int j = 1; j <= elnfa; j++)
	  nnums.Elem(elnv+elned+j) = nv+ned+fanums.Elem(j);
	nnums.Elem(elnv+elned+elnfa+1) = nv+ned+nfa+i;

	for (int j = 0; j < nnums.Size(); j++)
	  cluster_reps.Elem(nnums[j]) = nnums[j];
	  }
	*/
	netgen.GlobalMembers.ParallelForRange(new TaskManager(tm), ne, (uint begin, uint end) =>
	{
		 Array<int> nnums = new Array<int>();
		 Array<int> ednums = new Array<int>();
		 Array<int> fanums = new Array<int>();
		 for (int i = (int)(begin + 1); i <= end; i++)
		 {
			 Element el = mesh.VolumeElement(i);
			 ELEMENT_TYPE typ = el.GetType();

			 top.GetElementEdges(i, ednums);
			 top.GetElementFaces(i, fanums);

			 int elnv = top.GetNVertices(typ);
			 int elned = ednums.Size();
			 int elnfa = fanums.Size();

			 nnums.SetSize(elnv + elned + elnfa + 1);
			 for (int j = 1; j <= elnv; j++)
			 {
			   nnums.Elem(j) = el.PNum(j) + 1 - PointIndex.BASE;
			 }
			 for (int j = 1; j <= elned; j++)
			 {
			   nnums.Elem(elnv + j) = nv + ednums.Elem(j);
			 }
			 for (int j = 1; j <= elnfa; j++)
			 {
			   nnums.Elem(elnv + elned + j) = nv + ned + fanums.Elem(j);
			 }
			 nnums.Elem(elnv + elned + elnfa + 1) = nv + ned + nfa + i;

			 for (int j = 0; j < nnums.Size(); j++)
			 {
			   cluster_reps.Elem(nnums[j]) = nnums[j];
			 }
		 }
	});

	NgProfiler.StopTimer(Update_timer1);
	NgProfiler.StartTimer(Update_timer2);
	/*
	for (int i = 1; i <= nse; i++)
	  {
	const Element2d & el = mesh.SurfaceElement(i);
	ELEMENT_TYPE typ = el.GetType();
	
	top.GetSurfaceElementEdges (i, ednums);
	int fanum = top.GetSurfaceElementFace (i);
	
	int elnv = top.GetNVertices (typ);
	int elned = ednums.Size();

	nnums.SetSize(elnv+elned+1);
	for (int j = 1; j <= elnv; j++)
	  nnums.Elem(j) = el.PNum(j)+1-PointIndex::BASE;
	for (int j = 1; j <= elned; j++)
	  nnums.Elem(elnv+j) = nv+ednums.Elem(j);
	nnums.Elem(elnv+elned+1) = fanum;

	for (int j = 0; j < nnums.Size(); j++)
	  cluster_reps.Elem(nnums[j]) = nnums[j];
	  }
	*/
	netgen.GlobalMembers.ParallelForRange(new TaskManager(tm), (uint)nse, (uint begin, uint end) =>
	{
		 ArrayMem<int,9> nnums = new ArrayMem<int,9>();
		 ArrayMem<int,9> ednums = new ArrayMem<int,9>();
		 for (int i = (int)(begin + 1); i <= end; i++)
		 {
			 Element2d el = mesh.SurfaceElement(i);
			 ELEMENT_TYPE typ = el.GetType();

			 top.GetSurfaceElementEdges(i, ednums);
			 int fanum = top.GetSurfaceElementFace(i);

			 int elnv = top.GetNVertices(typ);
			 int elned = ednums.Size();

			 nnums.SetSize(elnv + elned + 1);
			 for (int j = 1; j <= elnv; j++)
			 {
			   nnums.Elem(j) = el.PNum(j) + 1 - PointIndex.BASE;
			 }
			 for (int j = 1; j <= elned; j++)
			 {
			   nnums.Elem(elnv + j) = nv + ednums.Elem(j);
			 }
			 nnums.Elem(elnv + elned + 1) = fanum;

			 for (int j = 0; j < nnums.Size(); j++)
			 {
			   cluster_reps.Elem(nnums[j]) = nnums[j];
			 }
		 }
	});


	NgProfiler.StopTimer(Update_timer2);
	NgProfiler.StartTimer(Update_timer3);


	int[] hex_cluster = {1, 2, 3, 4, 1, 2, 3, 4, 5, 6, 7, 8, 5, 6, 7, 8, 1, 2, 3, 4, 9, 9, 5, 8, 6, 7, 9};

	int[] prism_cluster = {1, 2, 3, 1, 2, 3, 4, 5, 6, 4, 5, 6, 3, 1, 2, 7, 7, 4, 5, 6, 7};

	int[] pyramid_cluster = {1, 2, 2, 1, 3, 4, 2, 1, 4, 6, 5, 5, 6, 7, 5, 7, 6, 4, 7};
	int[] tet_cluster14 = {1, 2, 3, 1, 1, 4, 5, 4, 5, 6, 7, 5, 4, 7, 7};

	int[] tet_cluster12 = {1, 1, 2, 3, 4, 4, 5, 1, 6, 6, 7, 7, 4, 6, 7};

	int[] tet_cluster13 = {1, 2, 1, 3, 4, 6, 4, 5, 1, 5, 7, 4, 7, 5, 7};

	int[] tet_cluster23 = {2, 1, 1, 3, 6, 5, 5, 4, 4, 1, 5, 7, 7, 4, 7};

	int[] tet_cluster24 = {2, 1, 3, 1, 4, 1, 5, 4, 6, 5, 5, 7, 4, 7, 7};

	int[] tet_cluster34 = {2, 3, 1, 1, 4, 5, 1, 6, 4, 5, 5, 4, 7, 7, 7};

	int cnt = 0;

	do
	{
		tracer("update cluster, identify", false);
	cnt++;
	changed = 0;

	for (int i = 1; i <= ne; i++)
	{
		Element el = mesh.VolumeElement(i);
		ELEMENT_TYPE typ = el.GetType();

//C++ TO C# CONVERTER TODO TASK: C# does not have an equivalent to pointers to value types:
//ORIGINAL LINE: const int * clustertab = null;
		int clustertab = null;
		switch (typ)
		{
		  case PRISM:
		  case PRISM12:
		clustertab = prism_cluster;
		break;
		  case HEX:
		clustertab = hex_cluster;
		break;
		  case PYRAMID:
		clustertab = pyramid_cluster;
		break;
		  case TET:
		  case TET10:
		if (cluster_reps.Get(el.PNum(1) + 1 - PointIndex.BASE) == cluster_reps.Get(el.PNum(2) + 1 - PointIndex.BASE))
		{
		  clustertab = tet_cluster12;
		}
		else if (cluster_reps.Get(el.PNum(1) + 1 - PointIndex.BASE) == cluster_reps.Get(el.PNum(3) + 1 - PointIndex.BASE))
		{
		  clustertab = tet_cluster13;
		}
		else if (cluster_reps.Get(el.PNum(1) + 1 - PointIndex.BASE) == cluster_reps.Get(el.PNum(4) + 1 - PointIndex.BASE))
		{
		  clustertab = tet_cluster14;
		}
		else if (cluster_reps.Get(el.PNum(2) + 1 - PointIndex.BASE) == cluster_reps.Get(el.PNum(3) + 1 - PointIndex.BASE))
		{
		  clustertab = tet_cluster23;
		}
		else if (cluster_reps.Get(el.PNum(2) + 1 - PointIndex.BASE) == cluster_reps.Get(el.PNum(4) + 1 - PointIndex.BASE))
		{
		  clustertab = tet_cluster24;
		}
		else if (cluster_reps.Get(el.PNum(3) + 1 - PointIndex.BASE) == cluster_reps.Get(el.PNum(4) + 1 - PointIndex.BASE))
		{
		  clustertab = tet_cluster34;
		}

		else
		{
		  clustertab = null;
		}
		break;
		  default:
		clustertab = null;
		break;
		}

		if (clustertab != 0)
		{
				top.GetElementEdges(i, ednums);
				top.GetElementFaces(i, fanums);

				int elnv = top.GetNVertices(typ);
				int elned = ednums.Size();
				int elnfa = fanums.Size();

				nnums.SetSize(elnv + elned + elnfa + 1);
				for (int j = 1; j <= elnv; j++)
				{
				  nnums.Elem(j) = el.PNum(j) + 1 - PointIndex.BASE;
				}
				for (int j = 1; j <= elned; j++)
				{
				  nnums.Elem(elnv + j) = nv + ednums.Elem(j);
				}
				for (int j = 1; j <= elnfa; j++)
				{
				  nnums.Elem(elnv + elned + j) = nv + ned + fanums.Elem(j);
				}
				nnums.Elem(elnv + elned + elnfa + 1) = nv + ned + nfa + i;



		  for (int j = 0; j < nnums.Size(); j++)
		  {
		for (int k = 0; k < j; k++)
		{
		  if (clustertab[j] == clustertab[k])
		  {
			  int jj = nnums[j];
			  int kk = nnums[k];

			  if (cluster_reps.Get(kk) < cluster_reps.Get(jj))
			  {
			swap(jj,kk);
			  }

			  if (cluster_reps.Get(jj) < cluster_reps.Get(kk))
			  {
			  /*
			  cluster_reps.Elem(kk) = cluster_reps.Get(jj);
			  changed = 1;
			  */

			  int rep = cluster_reps.Get(jj);
			  int next = cluster_reps.Get(kk);
			  do
			  {
				  int cur = next;
				  next = llist.Elem(next);

				  cluster_reps.Elem(cur) = rep;
				  llist.Elem(cur) = llist.Elem(rep);
				  llist.Elem(rep) = cur;
			  } while (next != 0);
			  changed = 1;
			  }
		  }
		}
		  }
		}

		/*
		  if (clustertab)
		  {
		  if (typ == PYRAMID)
		  (*testout) << "pyramid";
		  else if (typ == PRISM || typ == PRISM12)
		  (*testout) << "prism";
		  else if (typ == TET || typ == TET10)
		  (*testout) << "tet";
		  else
		  (*testout) << "unknown type" << endl;

		  (*testout) << ", nnums  = ";
		  for (j = 0; j < nnums.Size(); j++)
		  (*testout) << "node " << j << " = " << nnums[j] << ", rep = "
		  << cluster_reps.Get(nnums[j]) << endl;
		  }
		*/
	}
		tracer("update cluster, identify", true);
	} while (changed != 0);
	NgProfiler.StopTimer(Update_timer3);
	/*
	  (*testout) << "cluster reps:" << endl;
	  for (i = 1; i <= cluster_reps.Size(); i++)
	  {
	  (*testout) << i << ": ";
	  if (i <= nv)
	  (*testout) << "v" << i << " ";
	  else if (i <= nv+ned)
	  (*testout) << "e" << i-nv << " ";
	  else if (i <= nv+ned+nfa)
	  (*testout) << "f" << i-nv-ned << " ";
	  else
	  (*testout) << "c" << i-nv-ned-nfa << " ";
	  (*testout) << cluster_reps.Get(i) << endl;
	  }
	*/
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetVertexRepresentant(int vnr) const
  public int GetVertexRepresentant(int vnr)
  {
	  return cluster_reps.Get(vnr);
  }
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetEdgeRepresentant(int ednr) const
  public int GetEdgeRepresentant(int ednr)
  {
	  return cluster_reps.Get(nv + ednr);
  }
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetFaceRepresentant(int fnr) const
  public int GetFaceRepresentant(int fnr)
  {
	  return cluster_reps.Get(nv + ned + fnr);
  }
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetElementRepresentant(int enr) const
  public int GetElementRepresentant(int enr)
  {
	  return cluster_reps.Get(nv + ned + nfa + enr);
  }
}



