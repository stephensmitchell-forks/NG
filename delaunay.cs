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


  public class DelaunayTet
  {
	private PointIndex[] pnums = Arrays.InitializeWithDefaultInstances<PointIndex>(4);
	private int[] nb = new int[4];

	public DelaunayTet()
	{
		;
	}

	public DelaunayTet(DelaunayTet el)
	{
	  for (int i = 0; i < 4; i++)
	  {
	pnums[i] = el[i];
	  }
	}

	public DelaunayTet(Element el)
	{
	  for (int i = 0; i < 4; i++)
	  {
	pnums[i] = el[i];
	  }
	}

	public PointIndex this [int i]
	{
		get
		{
			return new netgen.PointIndex(pnums[i]);
		}
		set
		{
			pnums[i] = value;
		}
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: PointIndex operator [] (int i) const
	public PointIndex this [int i]
	{
		get
		{
			return new netgen.PointIndex(pnums[i]);
		}
	}

	public int NB(int i)
	{
		return nb[i];
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int NB(int i) const
	public int NB(int i)
	{
		return nb[i];
	}


//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int FaceNr(INDEX_3 & face) const
	public int FaceNr(INDEX_3 face) // which face nr is it ?
	{
	  for (int i = 0; i < 3; i++)
	  {
	if (pnums[i] != face.I1() && pnums[i] != face.I2() && pnums[i] != face.I3())
	{
	  return i;
	}
	  }
	  return 3;
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: INDEX_3 GetFace(int i) const
	public INDEX_3 GetFace(int i)
	{
	  return new INDEX_3(pnums[GlobalMembers.deltetfaces[i][0]], pnums[GlobalMembers.deltetfaces[i][1]], pnums[GlobalMembers.deltetfaces[i][2]]);
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetFace(int i, Element2d & face) const
	public void GetFace(int i, Element2d face)
	{
	  // face.SetType(TRIG);
	  face[0] = pnums[GlobalMembers.deltetfaces[i][0]];
	  face[1] = pnums[GlobalMembers.deltetfaces[i][1]];
	  face[2] = pnums[GlobalMembers.deltetfaces[i][2]];
	}
  }







  /*
    Table to maintain neighbour elements
  */
  public class MeshNB
  {
	// face nodes -> one element
	private INDEX_3_CLOSED_HASHTABLE<int> faces = new INDEX_3_CLOSED_HASHTABLE<int>();

	// 
	private Array<DelaunayTet> tets;


	// estimated number of points
	public MeshNB(Array<DelaunayTet> atets, int np)
	{
		this.faces = new netgen.INDEX_3_CLOSED_HASHTABLE<int>(200);
		this.tets = new Array<DelaunayTet>(atets);
		;
	}

	// add element with 4 nodes
	public void Add(int elnr)
	{
	  DelaunayTet el = tets.Elem(elnr);

	  for (int i = 0; i < 4; i++)
	  {
	  INDEX_3 i3 = INDEX_3.Sort(el.GetFace(i));

	  int posnr;

	  if (!faces.PositionCreate(i3, ref posnr))
	  {
		  // face already in use
		  int othertet = faces.GetData(posnr);

		  el.NB(i) = othertet;
		  if (othertet != 0)
		  {
		  int fnr = tets.Get(othertet).FaceNr(i3);
		  tets.Elem(othertet).NB(fnr) = elnr;
		  }
	  }
	  else
	  {
		  faces.SetData(posnr, elnr);
		  el.NB(i) = 0;
	  }
	  }
	}

	// delete element with 4 nodes
	public void Delete(int elnr)
	{
	  DelaunayTet el = tets.Elem(elnr);
	  for (int i = 0; i < 4; i++)
	  {
	faces.Set(el.GetFace(i).Sort(), el.NB(i));
	  }
	}

	// get neighbour of element elnr in direction fnr 
	public int GetNB(int elnr, int fnr)
	{
	  return tets.Get(elnr).NB(fnr);
	}

	//
	public void ResetFaceHT(int size)
	{
	  faces.SetSize(size);
	}
  }





  /*
    connected lists of cosphereical elements
  */
  public class SphereList
  {
	private Array<int> links = new Array<int>();
	public SphereList()
	{
		;
	}

	public void AddElement(int elnr)
	{
	  if (elnr > links.Size())
	  {
	links.Append(1);
	  }
	  links.Elem(elnr) = elnr;
	}

	public void DeleteElement(int elnr)
	{
	  links.Elem(elnr) = 0;
	}

	public void ConnectElement(int eli, int toi)
	{
	  links.Elem(eli) = links.Get(toi);
	  links.Elem(toi) = eli;
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: void GetList(int eli, Array<int> & linked) const
	public void GetList(int eli, Array<int> linked)
	{
	  linked.SetSize(0);
	  int pi = eli;

	  do
	  {
	  if (pi <= 0 || pi > links.Size())
	  {
		  cerr << "link, error " << "\n";
		  cerr << "pi = " << pi << " linked.s = " << linked.Size() << "\n";
		  Environment.Exit(1);
	  }
	  if (linked.Size() > links.Size())
	  {
		  cerr << "links have loop" << "\n";
		  Environment.Exit(1);
	  }

	  linked.Append(pi);
	  pi = links.Get(pi);
	  } while (pi != eli);
	}
  }






//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
//C++ TO C# CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
  Timer Delaunay_t("Meshing3::Delaunay");

}
