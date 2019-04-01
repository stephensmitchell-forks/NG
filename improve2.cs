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



///
public class MeshOptimize2d
{
  private int faceindex;
  private int improveedges;
  private double metricweight;
  private int writestatus;

  ///
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  MeshOptimize2d();
  public virtual void Dispose()
  {
	  ;
  }
  ///
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  void ImproveMesh(Mesh mesh2d, MeshingParameters mp);
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  void ImproveMeshJacobian(Mesh mesh2d, MeshingParameters mp);
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  void ImproveVolumeMesh(Mesh mesh);
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  void ProjectBoundaryPoints(Array<int> surfaceindex, Array<Point<3> > from, Array<Point<3> > dest);

//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private int EdgeSwapping_timer = NgProfiler.CreateTimer("EdgeSwapping 2D");
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private int EdgeSwapping_timerstart = NgProfiler.CreateTimer("EdgeSwapping 2D start");

  public void EdgeSwapping(Mesh mesh, int usemetric)
  {
	if (!faceindex)
	{
	if (usemetric != 0)
	{
	  PrintMessage(3, "Edgeswapping, metric");
	}
	else
	{
	  PrintMessage(3, "Edgeswapping, topological");
	}

	for (faceindex = 1; faceindex <= mesh.GetNFD(); faceindex++)
	{
		EdgeSwapping(mesh, usemetric);

		if (multithread.terminate)
		{
		  throw new Exception("Meshing stopped");
		}
	}

	faceindex = 0;
	mesh.CalcSurfacesOfNode();
	return;
	}


//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static int timer = NgProfiler::CreateTimer("EdgeSwapping 2D");
	NgProfiler.RegionTimer reg1 = new NgProfiler.RegionTimer(EdgeSwapping_timer);

//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static int timerstart = NgProfiler::CreateTimer("EdgeSwapping 2D start");
	NgProfiler.StartTimer(EdgeSwapping_timerstart);


	Array<SurfaceElementIndex> seia = new Array<SurfaceElementIndex>();
	mesh.GetSurfaceElementsOfFace(faceindex, seia);

	for (int i = 0; i < seia.Size(); i++)
	{
	  if (mesh[seia[i]].GetNP() != 3)
	  {
	  GenericImprove(mesh);
	  return;
	  }
	}

	int surfnr = mesh.GetFaceDescriptor(faceindex).SurfNr();

	Array<Neighbour> neighbors = new Array<Neighbour>(mesh.GetNSE());
	INDEX_2_HASHTABLE<trionedge> other = new INDEX_2_HASHTABLE<trionedge>(seia.Size() + 2);


	Array<char> swapped = new Array<char>(mesh.GetNSE());
	Array<int,PointIndex.BASE> pdef = new Array<int,PointIndex.BASE>(mesh.GetNP());
	Array<double,PointIndex.BASE> pangle = new Array<double,PointIndex.BASE>(mesh.GetNP());


	// int e;
	// double d;
	// Vec3d nv1, nv2;

	// double loch(-1);
	double[] minangle = {0, 1.481, 2.565, 3.627, 4.683, 5.736, 7, 9};


	for (int i = 0; i < seia.Size(); i++)
	{
	Element2d sel = mesh[seia[i]];
	for (int j = 0; j < 3; j++)
	{
	  pangle[sel[j]] = 0.0;
	}
	}
	// pangle = 0;

	for (int i = 0; i < seia.Size(); i++)
	{
	Element2d sel = mesh[seia[i]];
	for (int j = 0; j < 3; j++)
	{
		POINTTYPE typ = mesh[sel[j]].Type();
		if (typ == FIXEDPOINT || typ == EDGEPOINT)
		{
		pangle[sel[j]] += Angle(mesh[sel[(j + 1) % 3]] - mesh[sel[j]], mesh[sel[(j + 2) % 3]] - mesh[sel[j]]);
		}
	}
	}

	// for (PointIndex pi = PointIndex::BASE;
	// pi < mesh.GetNP()+PointIndex::BASE; pi++)

	// pdef = 0;
	for (int i = 0; i < seia.Size(); i++)
	{
	Element2d sel = mesh[seia[i]];
	for (int j = 0; j < 3; j++)
	{
		PointIndex pi = sel[j];
		if (mesh[pi].Type() == INNERPOINT || mesh[pi].Type() == SURFACEPOINT)
		{
		  pdef[pi] = -6;
		}
		else
		{
		  for (int j = 0; j < 8; j++)
		  {
		if (pangle[pi] >= minangle[j])
		{
		  pdef[pi] = -1 - j;
		}
		  }
		}
	}
	}

	for (int i = 0; i < seia.Size(); i++)
	{
	Element2d sel = mesh[seia[i]];
	for (int j = 0; j < 3; j++)
	{
	  pdef[sel[j]]++;
	}
	}

	for (int i = 0; i < seia.Size(); i++)
	{
	for (int j = 0; j < 3; j++)
	{
		neighbors[seia[i]].SetNr(j, -1);
		neighbors[seia[i]].SetOrientation(j, 0);
	}
	}

	/*
	  Array<Vec3d> normals(mesh.GetNP());
	  for (i = 1; i <= mesh.GetNSE(); i++)
	  {
	  Element2d & hel = mesh.SurfaceElement(i);
	  if (hel.GetIndex() == faceindex)
	  for (k = 1; k <= 3; k++)
	  {
	  int pi = hel.PNum(k);
	  SelectSurfaceOfPoint (mesh.Point(pi), hel.GeomInfoPi(k));
	  int surfi = mesh.GetFaceDescriptor(faceindex).SurfNr();
	  GetNormalVector (surfi, mesh.Point(pi), normals.Elem(pi));
	  normals.Elem(pi) /= normals.Elem(pi).Length();
	  }
	  }
	*/

	for (int i = 0; i < seia.Size(); i++)
	{
	Element2d sel = mesh[seia[i]];

	for (int j = 0; j < 3; j++)
	{
		PointIndex pi1 = sel.PNumMod(j + 2);
		PointIndex pi2 = sel.PNumMod(j + 3);

		//	    double loch = mesh.GetH(mesh[pi1]);

		INDEX_2 edge = new INDEX_2(pi1, pi2);
		edge.Sort();

		if (mesh.IsSegment(pi1, pi2))
		{
		  continue;
		}

		/*
		  if (segments.Used (edge))
		  continue;
		*/
		INDEX_2 ii2 = new INDEX_2(pi1, pi2);
		if (other.Used(ii2))
		{
		// INDEX_2 i2s(ii2);
		// i2s.Sort();

		int i2 = other.Get(ii2).tnr;
		int j2 = other.Get(ii2).sidenr;

		neighbors[seia[i]].SetNr(j, i2);
		neighbors[seia[i]].SetOrientation(j, j2);
		neighbors[i2].SetNr(j2, seia[i]);
		neighbors[i2].SetOrientation(j2, j);
		}
		else
		{
		other.Set(new INDEX_2(pi2, pi1), new trionedge(seia[i], j));
		}
	}
	}

	for (int i = 0; i < seia.Size(); i++)
	{
	  swapped[seia[i]] = 0;
	}

	NgProfiler.StopTimer(EdgeSwapping_timerstart);



	int t = 4;
	int done = 0;
	while (done == 0 && t >= 2)
	{
	for (int i = 0; i < seia.Size(); i++)
	{
		SurfaceElementIndex t1 = seia[i];

		if (mesh[t1].IsDeleted())
		{
		  continue;
		}

		if (mesh[t1].GetIndex() != faceindex)
		{
		  continue;
		}

		if (multithread.terminate)
		{
		  throw new Exception("Meshing stopped");
		}

		for (int o1 = 0; o1 < 3; o1++)
		{
		bool should;


		SurfaceElementIndex t2 = neighbors[t1].GetNr(o1);
		int o2 = neighbors[t1].GetOrientation(o1);

		if (t2 == -1)
		{
			continue;
		}
		if (swapped[t1] || swapped[t2])
		{
			continue;
		}


		PointIndex pi1 = mesh[t1].PNumMod(o1 + 1 + 1);
		PointIndex pi2 = mesh[t1].PNumMod(o1 + 1 + 2);
		PointIndex pi3 = mesh[t1].PNumMod(o1 + 1);
		PointIndex pi4 = mesh[t2].PNumMod(o2 + 1);

		PointGeomInfo gi1 = mesh[t1].GeomInfoPiMod(o1 + 1 + 1);
		PointGeomInfo gi2 = mesh[t1].GeomInfoPiMod(o1 + 1 + 2);
		PointGeomInfo gi3 = mesh[t1].GeomInfoPiMod(o1 + 1);
		PointGeomInfo gi4 = mesh[t2].GeomInfoPiMod(o2 + 1);

		bool allowswap = true;

		Vec < 3> auxvec1 = mesh[pi3] - mesh[pi4];
		Vec < 3> auxvec2 = mesh[pi1] - mesh[pi4];

		allowswap = allowswap && ngsimd.GlobalMembers.fabs(1.0 - (auxvec1 * auxvec2) / (auxvec1.Length() * auxvec2.Length())) > 1e-4;

		if (!allowswap)
		{
		  continue;
		}

		// normal of new
		Vec < 3> nv1 = netgen.GlobalMembers.Cross(auxvec1, auxvec2);

		auxvec1 = new mesh.Point(pi4) - new mesh.Point(pi3);
		auxvec2 = new mesh.Point(pi2) - new mesh.Point(pi3);
		allowswap = allowswap && ngsimd.GlobalMembers.fabs(1.0 - (auxvec1 * auxvec2) / (auxvec1.Length() * auxvec2.Length())) > 1e-4;


		if (!allowswap)
		{
		  continue;
		}

		Vec < 3> nv2 = netgen.GlobalMembers.Cross(auxvec1, auxvec2);


		// normals of original
		Vec < 3> nv3 = netgen.GlobalMembers.Cross(mesh[pi1] - mesh[pi4], mesh[pi2] - mesh[pi4]);
		Vec < 3> nv4 = netgen.GlobalMembers.Cross(mesh[pi2] - mesh[pi3], mesh[pi1] - mesh[pi3]);

		nv3 *= -1;
		nv4 *= -1;
		nv3.Normalize();
		nv4.Normalize();

		nv1.Normalize();
		nv2.Normalize();

		Vec < 3> nvp3, nvp4;
		SelectSurfaceOfPoint(new mesh.Point(pi3), gi3);
		GetNormalVector(surfnr, new mesh.Point(pi3), gi3, nvp3);

		nvp3.Normalize();

		SelectSurfaceOfPoint(new mesh.Point(pi4), gi4);
		GetNormalVector(surfnr, new mesh.Point(pi4), gi4, nvp4);

		nvp4.Normalize();



		double critval = ngsimd.GlobalMembers.cos(DefineConstants.M_PI / 6); // 30 degree
		allowswap = allowswap && (nv1 * nvp3 > critval) && (nv1 * nvp4 > critval) && (nv2 * nvp3 > critval) && (nv2 * nvp4 > critval) && (nvp3 * nv3 > critval) && (nvp4 * nv4 > critval);


		double horder = netgen.GlobalMembers.Dist(new mesh.Point(pi1), new mesh.Point(pi2));

		if (nv1.Length() > 1e-3 * horder * horder && nv2.Length() > 1e-3 * horder * horder && allowswap)
		{
			if (usemetric == 0)
			{
			int e = pdef[pi1] + pdef[pi2] - pdef[pi3] - pdef[pi4];
			double d = netgen.GlobalMembers.Dist2(new mesh.Point(pi1), new mesh.Point(pi2)) - netgen.GlobalMembers.Dist2(new mesh.Point(pi3), new mesh.Point(pi4));

			should = e >= t && (e > 2 || d > 0);
			}
			else
			{
			double loch = mesh.GetH(mesh[pi1]);
			should = CalcTriangleBadness(new mesh.Point(pi4), new mesh.Point(pi3), new mesh.Point(pi1), metricweight, loch) + CalcTriangleBadness(new mesh.Point(pi3), new mesh.Point(pi4), new mesh.Point(pi2), metricweight, loch) < CalcTriangleBadness(new mesh.Point(pi1), new mesh.Point(pi2), new mesh.Point(pi3), metricweight, loch) + CalcTriangleBadness(new mesh.Point(pi2), new mesh.Point(pi1), new mesh.Point(pi4), metricweight, loch);

			}

			if (allowswap)
			{
			Element2d sw1 = new Element2d(pi4, pi3, pi1);
			Element2d sw2 = new Element2d(pi3, pi4, pi2);

			int legal1 = mesh.LegalTrig(mesh.SurfaceElement(t1)) + mesh.LegalTrig(mesh.SurfaceElement(t2));
			int legal2 = mesh.LegalTrig(sw1) + mesh.LegalTrig(sw2);

			if (legal1 < legal2)
			{
				should = true;
			}
			if (legal2 < legal1)
			{
				should = false;
			}
			}

			if (should)
			{
			// do swapping !

			done = 1;

			mesh[t1].PNum(1) = pi1;
			mesh[t1].PNum(2) = pi4;
			mesh[t1].PNum(3) = pi3;

			mesh[t2].PNum(1) = pi2;
			mesh[t2].PNum(2) = pi3;
			mesh[t2].PNum(3) = pi4;

			mesh[t1].GeomInfoPi(1) = gi1;
			mesh[t1].GeomInfoPi(2) = gi4;
			mesh[t1].GeomInfoPi(3) = gi3;

			mesh[t2].GeomInfoPi(1) = gi2;
			mesh[t2].GeomInfoPi(2) = gi3;
			mesh[t2].GeomInfoPi(3) = gi4;

			pdef[pi1]--;
			pdef[pi2]--;
			pdef[pi3]++;
			pdef[pi4]++;

			swapped[t1] = 1;
			swapped[t2] = 1;
			}
		}
		}
	}
	t--;
	}

	mesh.SetNextTimeStamp();
  }

//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private int CombineImprove_timer = NgProfiler.CreateTimer("Combineimprove 2D");
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private int CombineImprove_timerstart = NgProfiler.CreateTimer("Combineimprove 2D start");
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private int CombineImprove_timerstart1 = NgProfiler.CreateTimer("Combineimprove 2D start1");

  public void CombineImprove(Mesh mesh)
  {
	if (!faceindex)
	{
	PrintMessage(3, "Combine improve");

	for (faceindex = 1; faceindex <= mesh.GetNFD(); faceindex++)
	{
		CombineImprove(mesh);

		if (multithread.terminate)
		{
		  throw new Exception("Meshing stopped");
		}
	}
	faceindex = 0;
	return;
	}


//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static int timer = NgProfiler::CreateTimer("Combineimprove 2D");
	NgProfiler.RegionTimer reg = new NgProfiler.RegionTimer(CombineImprove_timer);

//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static int timerstart = NgProfiler::CreateTimer("Combineimprove 2D start");
	NgProfiler.StartTimer(CombineImprove_timerstart);


//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static int timerstart1 = NgProfiler::CreateTimer("Combineimprove 2D start1");
	NgProfiler.StartTimer(CombineImprove_timerstart1);



	// int i, j, k, l;
	// PointIndex pi;
	// SurfaceElementIndex sei;


	Array<SurfaceElementIndex> seia = new Array<SurfaceElementIndex>();
	mesh.GetSurfaceElementsOfFace(faceindex, seia);


	for (int i = 0; i < seia.Size(); i++)
	{
	  if (mesh[seia[i]].GetNP() != 3)
	  {
	return;
	  }
	}



	int surfnr = 0;
	if (faceindex)
	{
	  surfnr = mesh.GetFaceDescriptor(faceindex).SurfNr();
	}


	// PointIndex pi1, pi2;
	// MeshPoint p1, p2, pnew;
	double bad1;
	double bad2;
	Vec < 3> nv;

	int np = mesh.GetNP();
	//int nse = mesh.GetNSE();

	TABLE<SurfaceElementIndex,PointIndex.BASE> elementsonnode = new TABLE<SurfaceElementIndex,PointIndex.BASE>(np);
	Array<SurfaceElementIndex> hasonepi = new Array<SurfaceElementIndex>();
	Array<SurfaceElementIndex> hasbothpi = new Array<SurfaceElementIndex>();

	for (int i = 0; i < seia.Size(); i++)
	{
	Element2d el = mesh[seia[i]];
	for (int j = 0; j < el.GetNP(); j++)
	{
	  elementsonnode.Add(el[j], seia[i]);
	}
	}

	Array<bool,PointIndex.BASE> @fixed = new Array<bool,PointIndex.BASE>(np);
	@fixed = false;

	NgProfiler.StopTimer(CombineImprove_timerstart1);

	/*
	for (SegmentIndex si = 0; si < mesh.GetNSeg(); si++)
	  {
	INDEX_2 i2(mesh[si][0], mesh[si][1]);
	fixed[i2.I1()] = true;
	fixed[i2.I2()] = true;
	  }
	*/

	for (int i = 0; i < seia.Size(); i++)
	{
	Element2d sel = mesh[seia[i]];
	for (int j = 0; j < sel.GetNP(); j++)
	{
		PointIndex pi1 = sel.PNumMod(j + 2);
		PointIndex pi2 = sel.PNumMod(j + 3);
		if (mesh.IsSegment(pi1, pi2))
		{
		@fixed[pi1] = true;
		@fixed[pi2] = true;
		}
	}
	}



	for (int i = 0; i < mesh.LockedPoints().Size(); i++)
	{
	  @fixed[mesh.LockedPoints()[i]] = true;
	}



	Array<Vec < 3>,PointIndex.BASE> normals = new Array<Vec < 3>,PointIndex.BASE>(np);

	for (PointIndex pi = mesh.Points().Begin(); pi < mesh.Points().End(); pi++)
	{
	if (elementsonnode[pi].Size())
	{
		Element2d hel = mesh[elementsonnode[pi][0]];
		for (int k = 0; k < 3; k++)
		{
		  if (hel[k] == pi)
		  {
		  SelectSurfaceOfPoint(mesh[pi], hel.GeomInfoPi(k + 1));
		  GetNormalVector(surfnr, mesh[pi], hel.GeomInfoPi(k + 1), normals[pi]);
		  break;
		  }
		}
	}
	}

	NgProfiler.StopTimer(CombineImprove_timerstart);

	for (int i = 0; i < seia.Size(); i++)
	{
	SurfaceElementIndex sei = seia[i];
	Element2d elem = mesh[sei];
	if (elem.IsDeleted())
	{
		continue;
	}

	for (int j = 0; j < 3; j++)
	{
		PointIndex pi1 = elem[j];
		PointIndex pi2 = elem[(j + 1) % 3];

		if (pi1 < PointIndex.BASE || pi2 < PointIndex.BASE)
		{
		  continue;
		}

		/*
		  INDEX_2 i2(pi1, pi2);
		  i2.Sort();
		  if (segmentht.Used(i2))
		  continue;
		*/

		bool debugflag = false;

		if (debugflag)
		{
		(*testout) << "Combineimprove, face = " << faceindex << "pi1 = " << pi1 << " pi2 = " << pi2 << "\n";
		}

		/*
		// save version:
		if (fixed.Get(pi1) || fixed.Get(pi2))
		continue;
		if (pi2 < pi1) swap (pi1, pi2);
		*/

		// more general
		if (@fixed[pi2])
		{
		  netgen.GlobalMembers.Swap(ref pi1, ref pi2);
		}

		if (@fixed[pi2])
		{
		  continue;
		}

		double loch = mesh.GetH(mesh[pi1]);

		INDEX_2 si2 = new INDEX_2(pi1, pi2);
		si2.Sort();

		/*
		  if (edgetested.Used (si2))
		  continue;
		  edgetested.Set (si2, 1);
		*/

		hasonepi.SetSize(0);
		hasbothpi.SetSize(0);

		for (int k = 0; k < elementsonnode[pi1].Size(); k++)
		{
		Element2d el2 = mesh[elementsonnode[pi1][k]];

		if (el2.IsDeleted())
		{
			continue;
		}

		if (el2[0] == pi2 || el2[1] == pi2 || el2[2] == pi2)
		{
			hasbothpi.Append(elementsonnode[pi1][k]);
			nv = netgen.GlobalMembers.Cross(new Vec3d(mesh[el2[0]], mesh[el2[1]]), new Vec3d(mesh[el2[0]], mesh[el2[2]]));
		}
		else
		{
			hasonepi.Append(elementsonnode[pi1][k]);
		}
		}


		Element2d hel = mesh[hasbothpi[0]];
		for (int k = 0; k < 3; k++)
		{
		  if (hel[k] == pi1)
		  {
		  SelectSurfaceOfPoint(mesh[pi1], hel.GeomInfoPi(k + 1));
		  GetNormalVector(surfnr, mesh[pi1], hel.GeomInfoPi(k + 1), nv);
		  break;
		  }
		}

		//	  nv = normals.Get(pi1);



		for (int k = 0; k < elementsonnode[pi2].Size(); k++)
		{
		Element2d el2 = mesh[elementsonnode[pi2][k]];
		if (el2.IsDeleted())
		{
			continue;
		}

		if (el2[0] == pi1 || el2[1] == pi1 || el2[2] == pi1)
		{
		  ;
		}
		else
		{
		  hasonepi.Append(elementsonnode[pi2][k]);
		}
		}

		bad1 = 0;
		int illegal1 = 0;
		int illegal2 = 0;
		for (int k = 0; k < hasonepi.Size(); k++)
		{
		Element2d el = mesh[hasonepi[k]];
		bad1 += CalcTriangleBadness(mesh[el[0]], mesh[el[1]], mesh[el[2]], nv, -1, loch);
		illegal1 += 1 - mesh.LegalTrig(el);
		}

		for (int k = 0; k < hasbothpi.Size(); k++)
		{
		Element2d el = mesh[hasbothpi[k]];
		bad1 += CalcTriangleBadness(mesh[el[0]], mesh[el[1]], mesh[el[2]], nv, -1, loch);
		illegal1 += 1 - mesh.LegalTrig(el);
		}
		bad1 /= (hasonepi.Size() + hasbothpi.Size());

		MeshPoint p1 = mesh[pi1];
		MeshPoint p2 = mesh[pi2];

		MeshPoint pnew = new MeshPoint(p1);
		mesh[pi1] = pnew;
		mesh[pi2] = pnew;

		bad2 = 0;
		for (int k = 0; k < hasonepi.Size(); k++)
		{
		Element2d el = mesh[hasonepi[k]];
		double err = CalcTriangleBadness(mesh[el[0]], mesh[el[1]], mesh[el[2]], nv, -1, loch);
		bad2 += err;

		Vec < 3> hnv = netgen.GlobalMembers.Cross(new Vec3d(mesh[el[0]], mesh[el[1]]), new Vec3d(mesh[el[0]], mesh[el[2]]));
		if (hnv * nv < 0)
		{
		  bad2 += 1e10;
		}

		for (int l = 0; l < 3; l++)
		{
		  if ((normals[el[l]] * nv) < 0.5)
		  {
			bad2 += 1e10;
		  }
		}

		illegal2 += 1 - mesh.LegalTrig(el);
		}
		bad2 /= hasonepi.Size();

		mesh[pi1] = p1;
		mesh[pi2] = p2;


		if (debugflag)
		{
		(*testout) << "bad1 = " << bad1 << ", bad2 = " << bad2 << "\n";
		}


		bool should = (bad2 < bad1 && bad2 < 1e4);
		if (bad2 < 1e4)
		{
		if (illegal1 > illegal2)
		{
			should = true;
		}
		if (illegal2 > illegal1)
		{
			should = false;
		}
		}


		if (should)
		{
				/*
				(*testout) << "combine !" << endl;
				(*testout) << "bad1 = " << bad1 << ", bad2 = " << bad2 << endl;
				(*testout) << "illegal1 = " << illegal1 << ", illegal2 = " << illegal2 << endl;
				(*testout) << "loch = " << loch << endl;
				*/

		mesh[pi1] = pnew;
		PointGeomInfo gi = new PointGeomInfo();
		// bool gi_set(false);


		Element2d el1p = new Element2d(null);
		int l = 0;
		while (mesh[elementsonnode[pi1][l]].IsDeleted() && l < elementsonnode.EntrySize(pi1))
		{
			l++;
		}
		if (l < elementsonnode.EntrySize(pi1))
		{
		  el1p = mesh[elementsonnode[pi1][l]];
		}
		else
		{
		  cerr << "OOPS!" << "\n";
		}

		for (l = 0; l < el1p.GetNP(); l++)
		{
		  if (el1p[l] == pi1)
		  {
			  gi = el1p.GeomInfoPi(l + 1);
			  // gi_set = true;
		  }
		}

		// (*testout) << "Connect point " << pi2 << " to " << pi1 << "\n";
		for (int k = 0; k < elementsonnode[pi2].Size(); k++)
		{
			Element2d el = mesh[elementsonnode[pi2][k]];
			if (el.IsDeleted())
			{
				continue;
			}
			elementsonnode.Add(pi1, elementsonnode[pi2][k]);

			bool haspi1 = false;
			for (l = 0; l < el.GetNP(); l++)
			{
			  if (el[l] == pi1)
			  {
			haspi1 = true;
			  }
			}
			if (haspi1)
			{
				continue;
			}

			for (int l = 0; l < el.GetNP(); l++)
			{
			if (el[l] == pi2)
			{
				el[l] = pi1;
				el.GeomInfoPi(l + 1) = gi;
			}

			@fixed[el[l]] = true;
			}
		}

		/*
		  for (k = 0; k < hasbothpi.Size(); k++)
		  {
		  cout << mesh[hasbothpi[k]] << endl;
		  for (l = 0; l < 3; l++)
		  cout << mesh[mesh[hasbothpi[k]][l]] << " ";
		  cout << endl;
		  }
		*/

		for (int k = 0; k < hasbothpi.Size(); k++)
		{
			mesh[hasbothpi[k]].Delete();
			/*
			  for (l = 0; l < 4; l++)
			  mesh[hasbothpi[k]][l] = PointIndex::BASE-1;
			*/
		}

		}
	}
	}

	//  mesh.Compress();
	mesh.SetNextTimeStamp();
  }

//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  void GenericImprove(Mesh mesh);


  public void SetFaceIndex(int fi)
  {
	  faceindex = fi;
  }
  public void SetImproveEdges(int ie)
  {
	  improveedges = ie;
  }
  public void SetMetricWeight(double mw)
  {
	  metricweight = mw;
  }
  public void SetWriteStatus(int ws)
  {
	  writestatus = ws;
  }



  ///
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  virtual void SelectSurfaceOfPoint(Point<3> p, PointGeomInfo gi);
  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual void ProjectPoint(INDEX, Point<3> &) const
  public virtual void ProjectPoint(INDEX UnnamedParameter, ref Point < 3>)
  {
  }

  /// project point, use gi as initial value, and compute new gi
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual int ProjectPointGI(INDEX surfind, Point<3> & p, PointGeomInfo & gi) const
  public virtual int ProjectPointGI(INDEX surfind, Point < 3> p, PointGeomInfo gi)
  {
	  ProjectPoint(surfind, p);
	  return CalcPointGeomInfo(surfind, gi, p);
  }

  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual void ProjectPoint2(INDEX, INDEX, Point<3> &) const
  public virtual void ProjectPoint2(INDEX UnnamedParameter, INDEX UnnamedParameter2, ref Point < 3>)
  {
  }

  /// liefert zu einem 3d-Punkt die geominfo (Dreieck) und liefert 1, wenn erfolgreich, 
  /// 0, wenn nicht (Punkt ausserhalb von chart)
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual int CalcPointGeomInfo(PointGeomInfo& gi, const Point<3> &) const
  public virtual int CalcPointGeomInfo(PointGeomInfo gi, Point < 3>)
  {
		gi.trignum = 1;
		return 1;
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual int CalcPointGeomInfo(int, PointGeomInfo& gi, const Point<3> & p3) const
  public virtual int CalcPointGeomInfo(int UnnamedParameter, PointGeomInfo gi, Point < 3> p3)
  {
		return CalcPointGeomInfo(gi, p3);
  }

  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual void GetNormalVector(INDEX surfind, const Point<3> & p, PointGeomInfo & gi, Vec<3> & n) const;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  virtual void GetNormalVector(INDEX surfind, Point<3> p, PointGeomInfo gi, Vec<3> n);
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual void GetNormalVector(INDEX surfind, const Point<3> & p, Vec<3> & n) const;
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  virtual void GetNormalVector(INDEX surfind, Point<3> p, Vec<3> n);

  public void CheckMeshApproximation(Mesh mesh)
  {
	// Check angles between elements and normals at corners
	/*

	int i, j;
	int ne = mesh.GetNSE();
	int surfnr;

	Vec3d n, ng;
	Array<Vec3d> ngs(3);

	(*mycout) << "Check Surface Approximation" << endl;
	(*testout) << "Check Surface Approximation" << endl;

	for (i = 1; i <= ne; i++)
	{
	const Element2d & el = mesh.SurfaceElement(i);
	surfnr = mesh.GetFaceDescriptor (el.GetIndex()).SurfNr();
	Vec3d n = Cross (mesh.Point (el.PNum(1)) - mesh.Point (el.PNum(2)),
	mesh.Point (el.PNum(1)) - mesh.Point (el.PNum(3)));
	n /= n.Length();

	for (j = 1; j <= el.GetNP(); j++)
	{
	SelectSurfaceOfPoint (mesh.Point(el.PNum(j)), el.GeomInfoPi(j));
	GetNormalVector (surfnr, mesh.Point(el.PNum(j)), ng);
	ng /= ng.Length();
	ngs.Elem(j) = ng;

	double angle =  (180.0 / M_PI) * Angle (n, ng);
	if (angle > 60)
	{
	(*testout) << "el " << i << " node " << el.PNum(j)
	<< "has angle = " << angle << endl;
	}
	}

	for (j = 1; j <= 3; j++)
	{
	double angle =  (180.0 / M_PI) * Angle (ngs.Get(j), ngs.Get(j%3+1));
	if (angle > 60)
	{
	(*testout) << "el " << i << " node-node "
	<< ngs.Get(j) << " - " << ngs.Get(j%3+1)
	<< " has angle = " << angle << endl;
	}
	}
	}
	*/
  }


  ///
//C++ TO C# CONVERTER TODO TASK: C# has no concept of a 'friend' class:
//  friend class Opti2SurfaceMinFunction;
  ///
//C++ TO C# CONVERTER TODO TASK: C# has no concept of a 'friend' class:
//  friend class Opti2EdgeMinFunction;
  ///
//C++ TO C# CONVERTER TODO TASK: C# has no concept of a 'friend' function:
//ORIGINAL LINE: friend double Opti2FunctionValueGrad(const Vector & x, Vector & grad);
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  double Opti2FunctionValueGrad(Vector x, Vector grad);
  ///
//C++ TO C# CONVERTER TODO TASK: C# has no concept of a 'friend' function:
//ORIGINAL LINE: friend double Opti2EdgeFunctionValueGrad(const Vector & x, Vector & grad);
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  double Opti2EdgeFunctionValueGrad(Vector x, Vector grad);



}

namespace netgen
{

  public class Neighbour
  {
	private int[] nr = new int[3];
	private int[] orient = new int[3];

	public Neighbour()
	{
		;
	}

	public void SetNr(int side, int anr)
	{
		nr[side] = anr;
	}
	public int GetNr(int side)
	{
		return nr[side];
	}

	public void SetOrientation(int side, int aorient)
	{
		orient[side] = aorient;
	}
	public int GetOrientation(int side)
	{
		return orient[side];
	}



	/*
	  void SetNr1 (int side, int anr) { nr[side-1] = anr; }
	  int GetNr1 (int side) { return nr[side-1]; }

	  void SetOrientation1 (int side, int aorient) { orient[side-1] = aorient; }
	  int GetOrientation1 (int side) { return orient[side-1]; }
	*/
  }




  public class trionedge
  {
	public int tnr;
	public int sidenr;

	public trionedge()
	{
		tnr = 0;
		sidenr = 0;
	}
	public trionedge(int atnr, int asidenr)
	{
		tnr = atnr;
		sidenr = asidenr;
	}
  }
}
