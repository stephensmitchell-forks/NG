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
/* File:   basegeom.hpp                                                   */
/* Author: Joachim Schoeberl                                              */
/* Date:   23. Aug. 09                                                    */
/**************************************************************************/


//C++ TO C# CONVERTER NOTE: C# has no need of forward class declarations:
//struct Tcl_Interp;

namespace netgen
{

//C++ TO C# CONVERTER WARNING: The original type declaration contained unconverted modifiers:
//ORIGINAL LINE: class DLL_HEADER NetgenGeometry
  public class NetgenGeometry : System.IDisposable
  {
	public virtual void Dispose()
	{
		;
	}

	public virtual int GenerateMesh(ref Mesh mesh, MeshingParameters mparam)
	{
	  if (mesh == null)
	  {
		  return 1;
	  }

	  if (mparam.perfstepsstart <= (int)MESHING_STEP.MESHCONST_MESHVOLUME)
	  {
	  multithread.task = "Volume meshing";

	  MESHING3_RESULT res = MeshVolume(mparam, mesh);

	  if (res != MESHING3_RESULT.MESHING3_OK)
	  {
		  return 1;
	  }

	  if (multithread.terminate)
	  {
		  return 0;
	  }

	  RemoveIllegalElements(mesh);
	  if (multithread.terminate)
	  {
		  return 0;
	  }

	  MeshQuality3d(mesh);
	  }


	  if (multithread.terminate || mparam.perfstepsend <= (int)MESHING_STEP.MESHCONST_MESHVOLUME)
	  {
		return 0;
	  }


	  if (mparam.perfstepsstart <= (int)MESHING_STEP.MESHCONST_OPTVOLUME)
	  {
	  multithread.task = "Volume optimization";

	  OptimizeVolume(mparam, mesh);
	  if (multithread.terminate)
	  {
		  return 0;
	  }
	  }

	  return 0;
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual const Refinement & GetRefinement() const
	public virtual Refinement GetRefinement()
	{
	  return new Refinement();
	}

  public virtual void DoArchive(Archive UnnamedParameter)
  {
//C++ TO C# CONVERTER TODO TASK: There is no C# equivalent to the C++ 'typeid' operator:
	  throw NgException("DoArchive not implemented for " + Demangle(typeid(this).name()));
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual void Save(string filename) const
	public virtual void Save(string filename)
	{
	  throw new Exception("Cannot save geometry - no geometry available");
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual void SaveToMeshFile(ostream &) const
	public virtual void SaveToMeshFile(ostream UnnamedParameter)
	{
		;
	}
  }





//C++ TO C# CONVERTER WARNING: The original type declaration contained unconverted modifiers:
//ORIGINAL LINE: class DLL_HEADER GeometryRegister
  public abstract class GeometryRegister : System.IDisposable
  {
	//DLL_HEADER Array<GeometryRegister*> geometryregister;

	public virtual void Dispose()
	{
		;
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual NetgenGeometry * Load(string filename) const = 0;
	public abstract NetgenGeometry Load(string filename);
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual NetgenGeometry * LoadFromMeshFile(istream &) const
	public virtual NetgenGeometry LoadFromMeshFile(istream UnnamedParameter)
	{
		return null;
	}
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual class VisualScene * GetVisualScene(const NetgenGeometry *) const
	public virtual VisualScene GetVisualScene(NetgenGeometry UnnamedParameter)
	{
		return null;
	}
	public virtual void SetParameters(Tcl_Interp UnnamedParameter)
	{
		;
	}
  }

//C++ TO C# CONVERTER WARNING: The original type declaration contained unconverted modifiers:
//ORIGINAL LINE: class DLL_HEADER GeometryRegisterArray : public Array<GeometryRegister*>
  public class GeometryRegisterArray : Array<GeometryRegister*>, System.IDisposable
  {
	public virtual void Dispose()
	{
	  for (int i = 0; i < Size(); i++)
	  {
		if (this[i] != null)
		{
			this[i].Dispose();
		}
	  }
	}

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual NetgenGeometry *LoadFromMeshFile(istream & ist) const
	public virtual NetgenGeometry LoadFromMeshFile(istream ist)
	{
	  for (int i = 0; i < Size(); i++)
	  {
		  NetgenGeometry hgeom = this[i].LoadFromMeshFile(ist);
		  if (hgeom != null)
		  {
			return hgeom;
		  }
	  }
	  return null;
	}
  }
}





