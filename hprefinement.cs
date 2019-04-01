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
/* File:   hprefinement.hh                                                */
/* Author: Joachim Schoeberl                                              */
/* Date:   27. Oct. 2000                                                  */
/**************************************************************************/

/*
  HP Refinement
*/




public enum HPREF_ELEMENT_TYPE
{
  HP_NONE = 0,

  HP_SEGM = 1,
  HP_SEGM_SINGCORNERL,
  HP_SEGM_SINGCORNERR,
  HP_SEGM_SINGCORNERS,

  HP_TRIG = 10,
  HP_TRIG_SINGCORNER,
  HP_TRIG_SINGCORNER12,
  HP_TRIG_SINGCORNER123,
  HP_TRIG_SINGCORNER123_2D, // not rotational symmetric
  HP_TRIG_SINGEDGE = 20,
  HP_TRIG_SINGEDGECORNER1, // E = 100, V = 100
  HP_TRIG_SINGEDGECORNER2, // E = 100, V = 010
  HP_TRIG_SINGEDGECORNER12, // E = 100, V = 110
  HP_TRIG_SINGEDGECORNER3,
  HP_TRIG_SINGEDGECORNER13,
  HP_TRIG_SINGEDGECORNER23,
  HP_TRIG_SINGEDGECORNER123,
  HP_TRIG_SINGEDGES = 30,
  HP_TRIG_SINGEDGES2,
  HP_TRIG_SINGEDGES3,
  HP_TRIG_SINGEDGES23,
  HP_TRIG_3SINGEDGES = 40,

  HP_QUAD = 50,
  HP_QUAD_SINGCORNER,
  HP_DUMMY_QUAD_SINGCORNER,
  HP_QUAD_SINGEDGE,
  HP_QUAD_0E_2VA, // V = 1100
  HP_QUAD_0E_2VB, // V = 1010
  HP_QUAD_0E_3V,
  HP_QUAD_0E_4V,

  // one edge: marked edge is always edge from vertex 1 to vertex 2 (E = 1000)
  HP_QUAD_1E_1VA, // vertex on beginning of edge: V = 1000
  HP_QUAD_1E_1VB, // vertex on end of edge: V = 0100
  HP_QUAD_1E_1VC, // V = 0010
  HP_QUAD_1E_1VD, // V = 0001

  HP_QUAD_1E_2VA, // V = 1100
  HP_QUAD_1E_2VB, // V = 1010
  HP_QUAD_1E_2VC, // V = 1001
  HP_QUAD_1E_2VD, // V = 0110
  HP_QUAD_1E_2VE, // V = 0101
  HP_QUAD_1E_2VF, // V = 0011

  HP_QUAD_1E_3VA, // V = 1110
  HP_QUAD_1E_3VB, // V = 1101
  HP_QUAD_1E_3VC, // V = 1011
  HP_QUAD_1E_3VD, // V = 0111

  HP_QUAD_1E_4V, // V = 1111


  HP_QUAD_2E, // E = 1001, V = 1000
  HP_QUAD_2E_1VA, // E = 1001, V = 1100
  HP_QUAD_2E_1VB, // E = 1001, V = 1010
  HP_QUAD_2E_1VC, // E = 1001, V = 1001
  HP_QUAD_2E_2VA, // E = 1001, V = 1110
  HP_QUAD_2E_2VB, // E = 1001, V = 1101
  HP_QUAD_2E_2VC, // E = 1001, V = 1011
  HP_QUAD_2E_3V, // E = 1001, V = 1111

  HP_QUAD_2EB_0V, // E = 1010, V = 0000
  HP_QUAD_2EB_1VA, // E = 1010, V = 1000
  HP_QUAD_2EB_1VB, // E = 1010, V = 0100
  HP_QUAD_2EB_2VA, // E = 1010, V = 1100
  HP_QUAD_2EB_2VB, // E = 1010, V = 1010
  HP_QUAD_2EB_2VC, // E = 1010, V = 1001
  HP_QUAD_2EB_2VD, // E = 1010, V = 0101
  HP_QUAD_2EB_3VA, // E = 1010, V = 1110
  HP_QUAD_2EB_3VB, // E = 1010, V = 1101

  HP_QUAD_2EB_4V,


  HP_QUAD_3E, // E = 1101, V = 1100
  HP_QUAD_3E_3VA, // E = 1101, V = 1110
  HP_QUAD_3E_3VB, // E = 1101, V = 1101
  HP_QUAD_3E_4V, // E = 1101, V = 1111

  HP_QUAD_4E,


  HP_TET = 100, // no singular vertex/edge
  HP_TET_0E_1V, // V1
  HP_TET_0E_2V, // V1,2
  HP_TET_0E_3V, // V1,2,3
  HP_TET_0E_4V, // V1,2,3,4
  HP_TET_1E_0V = 200, // E1-2
  HP_TET_1E_1VA, // V1
  HP_TET_1E_1VB, // V3
  HP_TET_1E_2VA, // V1,2
  HP_TET_1E_2VB, // V1,3
  HP_TET_1E_2VC, // V1,4
  HP_TET_1E_2VD, // V3,4
  HP_TET_1E_3VA, // V1,2,3
  HP_TET_1E_3VB, // V1,3,4
  HP_TET_1E_4V, // V1,2,3,4


  // 2 connected edges, additionally marked Vs
  HP_TET_2EA_0V = 220, // E1-2, E1-3
  HP_TET_2EA_1VA, // V2
  HP_TET_2EA_1VB, // V3
  HP_TET_2EA_1VC, // V4
  HP_TET_2EA_2VA, // V2,3
  HP_TET_2EA_2VB, // V2,4
  HP_TET_2EA_2VC, // V3,4
  HP_TET_2EA_3V, // V2,3,4

  // 2 opposite edges
  HP_TET_2EB_0V = 230, // E1-2, E3-4
  HP_TET_2EB_1V, // V1
  HP_TET_2EB_2VA, // V1,2
  HP_TET_2EB_2VB, // V1,3
  HP_TET_2EB_2VC, // V1,4
  HP_TET_2EB_3V, // V1,2,3
  HP_TET_2EB_4V, // V1,2,3,4

  HP_TET_3EA_0V = 400, // E1-2, E1-3, E1-4, 3 edges connected
  HP_TET_3EA_1V, // V2
  HP_TET_3EA_2V, // V2,3
  HP_TET_3EA_3V, // V2,3,4

  HP_TET_3EB_0V = 420, // E1-2, E1-4, E2-3  3 edges chain
  HP_TET_3EB_1V,
  HP_TET_3EB_2V,
  HP_TET_3EC_0V = 430, // 3 edges chain, alter
  HP_TET_3EC_1V, // 3 edges chain, alter
  HP_TET_3EC_2V, // 3 edges chain, alter


  HP_TET_1F_0E_0V = 500, // 1 singular face
  HP_TET_1F_0E_1VA, // 1 sing vertex in face (V2)
  HP_TET_1F_0E_1VB, // 1 sing vertex not in face (V1)
  HP_TET_1F_1EA_0V, // 1 sing edge not in face
  HP_TET_1F_1EB_0V, // 1 sing edge in face
  HP_TET_2F_0E_0V = 600, // 2 singular faces

  HP_PRISM = 1000,
  HP_PRISM_SINGEDGE,
  HP_PRISM_SINGEDGE_V12,
  HP_PRISM_SINGEDGE_H1,
  HP_PRISM_SINGEDGE_H12,

  HP_PRISM_1FA_0E_0V, // 1 singular trig face
  HP_PRISM_2FA_0E_0V, // 2 singular trig faces
  HP_PRISM_1FB_0E_0V, // 1 singular quad face  1-2-4-5

  HP_PRISM_1FB_1EA_0V, // 1 singular quad face, edge is 1-2
  HP_PRISM_1FA_1E_0V,
  HP_PRISM_2FA_1E_0V,
  HP_PRISM_1FA_1FB_0E_0V,
  HP_PRISM_2FA_1FB_0E_0V,
  HP_PRISM_1FA_1FB_1EA_0V,
  HP_PRISM_1FA_1FB_1EB_0V,
  HP_PRISM_2FA_1FB_1EA_0V,
  HP_PRISM_1FB_1EC_0V,
  HP_PRISM_1FA_1FB_1EC_0V,
  HP_PRISM_2FA_1FB_1EC_0V,
  HP_PRISM_1FB_2EA_0V,
  HP_PRISM_1FA_1FB_2EA_0V,
  HP_PRISM_2FA_1FB_2EA_0V,
  HP_PRISM_1FB_2EB_0V,
  HP_PRISM_1FA_1FB_2EB_0V,
  HP_PRISM_1FA_1FB_2EC_0V,
  HP_PRISM_2FA_1FB_2EB_0V,
  HP_PRISM_1FB_3E_0V,
  HP_PRISM_1FA_1FB_3E_0V,
  HP_PRISM_2FA_1FB_3E_0V,
  HP_PRISM_2FB_0E_0V,
  HP_PRISM_1FA_2FB_0E_0V,
  HP_PRISM_2FA_2FB_0E_0V,
  HP_PRISM_2FB_1EC_0V,
  HP_PRISM_1FA_2FB_1EC_0V,
  HP_PRISM_1FA_2FB_1EB_0V,
  HP_PRISM_2FA_2FB_1EC_0V,
  HP_PRISM_2FB_3E_0V,
  HP_PRISM_1FA_2FB_3E_0V,
  HP_PRISM_2FA_2FB_3E_0V,
  HP_PRISM_1FA_2E_0V,
  HP_PRISM_2FA_2E_0V,
  HP_PRISM_3E_0V,
  HP_PRISM_1FA_3E_0V,
  HP_PRISM_2FA_3E_0V,
  HP_PRISM_3FB_0V,
  HP_PRISM_1FA_3FB_0V,
  HP_PRISM_2FA_3FB_0V,
  HP_PRISM_3E_4EH,



  /*  HP_PRISM_1FB_1EA_0V,     // 1 singular quad face, edge is 1-4
  HP_PRISM_1FB_1EB_0V,     // 1 singular quad face, edge is 2-5
  HP_PRISM_2F_0E_0V,      // 2 singular quad faces
  */

  HP_PYRAMID = 2000,
  HP_PYRAMID_0E_1V,
  HP_PYRAMID_EDGES,
  HP_PYRAMID_1FB_0E_1VA, // 1 trig face, top vertex

  HP_HEX = 3000,
  HP_HEX_0E_1V,
  HP_HEX_1E_1V,
  HP_HEX_1E_0V,
  HP_HEX_3E_0V,
  HP_HEX_1F_0E_0V,
  HP_HEX_1FA_1FB_0E_0V
}



public class HPRef_Struct
{
  public HPREF_ELEMENT_TYPE geom = new HPREF_ELEMENT_TYPE();
  public int[] splitedges = new int[3];
  public int[] splitfaces = new int[4];
  public int[] splitelements = new int[5];
  public HPREF_ELEMENT_TYPE[] neweltypes;
  public int[] newels = new int[8];
}




public class HPRefElement
{
  private void Reset()
  {
	np = 8;
	for (int i = 0; i < 8; i++)
	{
	pnums[i] = -1;
	param[i][0] = param[i][1] = param[i][2] = 0;
	}
	domin = -1;
	domout = -1; // he:
	levelx = 0;
	levely = 0;
	levelz = 0;
  }

  public HPRefElement()
  {
	Reset();
  }

  public HPRefElement(Element el)
  {
	  this.np = el.GetNV();
	  this.index = el.GetIndex();
	  this.levelx = 0;
	  this.levely = 0;
	  this.levelz = 0;
	  this.type = HP_NONE;
	  this.domin = -1;
	  this.domout = -1;
	//Reset();
	for (int i = 0; i < np ; i++)
	{
	  pnums[i] = el[i];
	}

	Point3d[] points = MeshTopology.GetVertices(el.GetType());
	for (int i = 0;i < np;i++)
	{
	  for (int l = 0;l < 3;l++)
	  {
	param[i][l] = points[i].X(l + 1);
	  }
	}
  }

  public HPRefElement(Element2d el)
  {
	  this.levelx = 0;
	  this.levely = 0;
	  this.levelz = 0;
	  this.type = HP_NONE;
	  this.index = el.GetIndex();
	  this.np = el.GetNV();
	  this.domin = -1;
	  this.domout = -1;
	//Reset();

	for (int i = 0; i < np ; i++)
	{
	  pnums[i] = el[i];
	}

	Point3d[] points = MeshTopology.GetVertices(el.GetType());
	for (int i = 0;i < np;i++)
	{
	  for (int l = 0;l < 3;l++)
	  {
	param[i][l] = points[i].X(l + 1);
	  }
	}
  }

  public HPRefElement(Segment el)
  {
	  this.levelx = 0;
	  this.levely = 0;
	  this.levelz = 0;
	  this.type = HP_NONE;
	  this.np = 2;
	  this.domin = el.domin;
	  this.domout = el.domout;
	  this.singedge_left = el.singedge_left;
	  this.singedge_right = el.singedge_right;
	//Reset();
	for (int i = 0; i < np ; i++)
	{
	  pnums[i] = el[i];
	}
	Point3d[] points = MeshTopology.GetVertices(SEGMENT);
	for (int i = 0;i < np;i++)
	{
	  for (int l = 0;l < 3;l++)
	  {
		param[i][l] = points[i].X(l + 1);
	  }
	}

	/*
	for (int i=0; i<np; i++)
	{
	  param[i][0] = i;
	  param[i][1] = -1; param[i][2] = -1;
	}
	*/
  }

  public HPRefElement(HPRefElement el)

  {
	  this.np = el.np;
	  this.levelx = el.levelx;
	  this.levely = el.levely;
	  this.levelz = el.levelz;
	  this.type = el.type;
	  this.domin = el.domin;
	  this.domout = el.domout;
	  this.index = el.index;
	  this.coarse_elnr = el.coarse_elnr;
	  this.singedge_left = el.singedge_left;
	  this.singedge_right = el.singedge_right;
	//Reset();
	for (int i = 0; i < np ; i++)
	{
	pnums[i] = el[i];
	for (int l = 0; l < 3; l++)
	{
		param[i][l] = el.param[i][l];
	}
	}
  }

  public void SetType(HPREF_ELEMENT_TYPE t)
  {
	type = t;
	switch (type)
	{
	  case HP_SEGM:
		  np = 2;
		  break;
	  case HP_TRIG:
		  np = 3;
		  break;
	  case HP_QUAD:
		  np = 4;
		  break;
	  case HP_TET:
		  np = 4;
		  break;
	  case HP_PRISM:
		  np = 6;
		  break;
	  case HP_PYRAMID:
		  np = 5;
		  break;
	  case HP_HEX:
		  np = 8;
		  break;

	  default:
		cerr << "HPRefElement: illegal type " << type << "\n";
		throw new Exception("HPRefElement::SetType: illegal type");
	}

	for (int k = 0; k < 8;k++)
	{
	pnums[k] = 0;
	for (int l = 0; l < 3; l++)
	{
		  param[k][l] = 0.0;
	}
	}
  }

  // HPRefElement(HPRefElement & el, HPREF_ELEMENT_TYPE t); 

  /* HPRefElement(HPRefElement & el, HPREF_ELEMENT_TYPE t)
  { 
    type = t; 
    HPRef_Struct * hprs = Get_HPRef_Struct(t);
    for (int i=0; i<np ; i++) 
      {
	pnums[i] = el[i];
	for(int l=0; l<np; l++) param[i][l] = el.param[i][l]; 
      }
    switch(hprs->geom)
      {
      case HP_SEGM: np=2; sing_edge_left=0; sing_edge_right=0; break; 
      case HP_QUAD: np=4; break; 
      case HP_TRIG: np=3; break; 
      case HP_HEX: np=8; break; 
      case HP_PRISM: np=6; break;
      case HP_TET: np=4; break; 
      case HP_PYRAMID: np=5; break; 
      }
    index = el.index; 
    levelx = el.levelx; 
    levely = el.levely; 
    levelz = el.levelz; 
    type = el.type; 
    coarse_elnr = el.coarse_elnr;
    singedge_left = el.singedge_left; 
    singedge_right = el.singedge_left; 
    } */ 

  public HPREF_ELEMENT_TYPE type = new HPREF_ELEMENT_TYPE();
  public PointIndex[] pnums = Arrays.InitializeWithDefaultInstances<PointIndex>(8);
  public double[][] param = RectangularArrays.RectangularDoubleArray(8, 3);
  public int index;
  public int levelx;
  public int levely;
  public int levelz;
  public int np;
  public int coarse_elnr;
  public int domin; // he: needed for segment!! in 3d there should be surf1, surf2!!
  public int domout;
  // int coarse_hpelnr; 
  public PointIndex this[int i]
  {
	  get
	  {
		  return (pnums[i]);
	  }
	  set
	  {
		  (pnums[i]) = value;
	  }
  }
  public PointIndex PNumMod(int i)
  {
	  return pnums[(i - 1) % np];
  }
  public PointIndex PNum(int i)
  {
	  return pnums[(i - 1)];
  }
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: int GetIndex() const
  public int GetIndex()
  {
	  return index;
  }
  public double singedge_left;
  public double singedge_right;


  //  EdgePointGeomInfo epgeominfo[2];

}






