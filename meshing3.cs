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




public enum MESHING3_RESULT
{
  MESHING3_OK = 0,
  MESHING3_GIVEUP = 1,
  MESHING3_NEGVOL = 2,
  MESHING3_OUTERSTEPSEXCEEDED = 3,
  MESHING3_TERMINATE = 4,
  MESHING3_BADSURFACEMESH = 5
}


/// 3d volume mesh generation
public class Meshing3
{
  /// current state of front
  private AdFront3 adfront;
  /// 3d generation rules
  private Array<vnetrule> rules = new Array<vnetrule>();
  /// counts how often a rule is used
  private Array<int> ruleused = new Array<int>();
  private Array<int> canuse = new Array<int>();
  private Array<int> foundmap = new Array<int>();
  /// describes, why a rule is not applied
  private Array<char> problems = new Array<char>();
  /// tolerance criterion
  private double tolfak;
  /// 
  public Meshing3(string rulefilename)
  {
	tolfak = 1;

	LoadRules(rulefilename, null);
	adfront = new AdFront3();

	problems.SetSize(rules.Size());
	foundmap.SetSize(rules.Size());
	canuse.SetSize(rules.Size());
	ruleused.SetSize(rules.Size());

	for (int i = 1; i <= rules.Size(); i++)
	{
		problems.Elem(i) = new char[255];
		foundmap.Elem(i) = 0;
		canuse.Elem(i) = 0;
		ruleused.Elem(i) = 0;
	}
  }

  /// 
  public Meshing3(string[] rulep)
  {
	tolfak = 1;

	LoadRules(null, rulep);
	adfront = new AdFront3();

	problems.SetSize(rules.Size());
	foundmap.SetSize(rules.Size());
	canuse.SetSize(rules.Size());
	ruleused.SetSize(rules.Size());

	for (int i = 0; i < rules.Size(); i++)
	{
		problems[i] = new char[255];
		foundmap[i] = 0;
		canuse[i] = 0;
		ruleused[i] = 0;
	}
  }

  ///
  public virtual void Dispose()
  {
	if (adfront != null)
	{
		adfront.Dispose();
	}
	for (int i = 0; i < rules.Size(); i++)
	{
		Arrays.DeleteArray(problems[i]);
		if (rules[i] != null)
		{
			rules[i].Dispose();
		}
	}
  }

  ///
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  void LoadRules(string filename, string[] prules);
  ///
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
//C++ TO C# CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
  Timer GenerateMesh_t("Meshing3::GenerateMesh");
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
//C++ TO C# CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
  Timer GenerateMesh_meshing3_timer_a("Meshing3::GenerateMesh a", 2);
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
//C++ TO C# CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
  Timer GenerateMesh_meshing3_timer_b("Meshing3::GenerateMesh b", 2);
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
//C++ TO C# CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
  Timer GenerateMesh_meshing3_timer_c("Meshing3::GenerateMesh c", 1);
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
//C++ TO C# CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
  Timer GenerateMesh_meshing3_timer_d("Meshing3::GenerateMesh d", 2);

  public MESHING3_RESULT GenerateMesh(Mesh mesh, MeshingParameters mp)
  {
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static Timer t("Meshing3::GenerateMesh");
	RegionTimer reg = new RegionTimer(GenerateMesh_t);
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static Timer meshing3_timer_a("Meshing3::GenerateMesh a", 2);
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static Timer meshing3_timer_b("Meshing3::GenerateMesh b", 2);
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static Timer meshing3_timer_c("Meshing3::GenerateMesh c", 1);
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static Timer meshing3_timer_d("Meshing3::GenerateMesh d", 2);
	// static int meshing3_timer = NgProfiler::CreateTimer ("Meshing3::GenerateMesh");
	// static int meshing3_timer_a = NgProfiler::CreateTimer ("Meshing3::GenerateMesh a");
	// static int meshing3_timer_b = NgProfiler::CreateTimer ("Meshing3::GenerateMesh b");
	// static int meshing3_timer_c = NgProfiler::CreateTimer ("Meshing3::GenerateMesh c");
	// static int meshing3_timer_d = NgProfiler::CreateTimer ("Meshing3::GenerateMesh d");
	// NgProfiler::RegionTimer reg (meshing3_timer);


	Array<Point3d, PointIndex.BASE> locpoints = new Array<Point3d, PointIndex.BASE>(); // local points
	Array<MiniElement2d> locfaces = new Array<MiniElement2d>(); // local faces
	Array<PointIndex, PointIndex.BASE> pindex = new Array<PointIndex, PointIndex.BASE>(); // mapping from local to front point numbering
	Array<int, PointIndex.BASE> allowpoint = new Array<int, PointIndex.BASE>(); // point is allowed ?
	Array<int> findex = new Array<int>(); // mapping from local to front face numbering
	//INDEX_2_HASHTABLE<int> connectedpairs(100);    // connecgted pairs for prism meshing

	Array<Point3d, PointIndex.BASE> plainpoints = new Array<Point3d, PointIndex.BASE>(); // points in reference coordinates
	Array<int> delpoints = new Array<int>(); // points and lines to be deleted
	Array<int> delfaces = new Array<int>();
	Array<Element> locelements = new Array<Element>(); // new generated elements

	int j;
	int oldnp;
	int oldnf;
	int found;
	referencetransform trans = new referencetransform();
	int rotind;
	Point3d inp = new Point3d();
	float err;

	int locfacesplit; //index for faces in outer area

	bool loktestmode = false;

	int uselocalh = mp.uselocalh;

	// int giveuptol = mp.giveuptol; //
	MeshingStat3d stat = new MeshingStat3d(); // statistics
	int plotstat_oldne = -1;


	// for star-shaped domain meshing
	Array<MeshPoint, PointIndex.BASE> grouppoints = new Array<MeshPoint, PointIndex.BASE>();
	Array<MiniElement2d> groupfaces = new Array<MiniElement2d>();
	Array<PointIndex, PointIndex.BASE> grouppindex = new Array<PointIndex, PointIndex.BASE>();
	Array<int> groupfindex = new Array<int>();


	float minerr;
	int hasfound;
	double tetvol;
	// int giveup = 0;


	Array<Point3d> tempnewpoints = new Array<Point3d>();
	Array<MiniElement2d> tempnewfaces = new Array<MiniElement2d>();
	Array<int> tempdelfaces = new Array<int>();
	Array<Element> templocelements = new Array<Element>();


	stat.h = mp.maxh;

	adfront.SetStartFront(mp.baseelnp);


	found = 0;
	stat.vol0 = adfront.Volume();
	tetvol = 0;

	stat.qualclass = 1;

	while (true)
	{
		if (multithread.terminate)
		{
	  throw new Exception("Meshing stopped");
		}

		// break if advancing front is empty
		if (!mp.baseelnp && adfront.Empty())
		{
	  break;
		}

		// break, if advancing front has no elements with
		// mp.baseelnp nodes
		if (mp.baseelnp && adfront.Empty(mp.baseelnp))
		{
	  break;
		}

		locpoints.SetSize(0);
		locfaces.SetSize(0);
		locelements.SetSize(0);
		pindex.SetSize(0);
		findex.SetSize(0);

		INDEX_2_HASHTABLE<int> connectedpairs = new INDEX_2_HASHTABLE<int>(100); // connected pairs for prism meshing

		// select base-element (will be locface[1])
		// and get local environment of radius (safety * h)


		int baseelem = adfront.SelectBaseElement();
		if (mp.baseelnp && adfront.GetFace(baseelem).GetNP() != mp.baseelnp)
		{
		adfront.IncrementClass(baseelem);
		continue;
		}

		MiniElement2d bel = adfront.GetFace(baseelem);
		const Point < 3> p1 = adfront.GetPoint(bel[0]);
		const Point < 3> p2 = adfront.GetPoint(bel[1]);
		const Point < 3> p3 = adfront.GetPoint(bel[2]);


		Point < 3> pmid = netgen.GlobalMembers.Center(p1, p2, p3);

		double his = (netgen.GlobalMembers.Dist(p1, p2) + netgen.GlobalMembers.Dist(p1, p3) + netgen.GlobalMembers.Dist(p2, p3)) / 3;
		double hshould = mesh.GetH(pmid);

		if (adfront.GetFace(baseelem).GetNP() == 4)
		{
	  hshould = netgen.GlobalMembers.max2(his, hshould);
		}

		double hmax = (his > hshould) ? his : hshould;

		// qualclass should come from baseelem !!!!!
		double hinner = hmax * (1 + stat.qualclass);
		double houter = hmax * (1 + 2 * stat.qualclass);

		GenerateMesh_meshing3_timer_a.Start();
		stat.qualclass = adfront.GetLocals(baseelem, locpoints, locfaces, pindex, findex, connectedpairs, houter, hinner, locfacesplit);
		GenerateMesh_meshing3_timer_a.Stop();

		// (*testout) << "locfaces = " << endl << locfaces << endl;


		//loktestmode = 1;
		testmode = loktestmode; //changed
		// loktestmode = testmode =  (adfront->GetFace (baseelem).GetNP() == 4) && (rules.Size() == 5);

		loktestmode = stat.qualclass > 5;


		if (loktestmode)
		{
		(*testout) << "baseel = " << baseelem << ", ind = " << findex.Get(1) << "\n";
			int pi1 = pindex[locfaces[0].PNum(1)];
			int pi2 = pindex[locfaces[0].PNum(2)];
			int pi3 = pindex[locfaces[0].PNum(3)];
		(*testout) << "pi = " << pi1 << ", " << pi2 << ", " << pi3 << "\n";
		}


		if (testmode)
		{
		(*testout) << "baseelem = " << baseelem << " qualclass = " << stat.qualclass << "\n";
		(*testout) << "locpoints = " << "\n" << locpoints << "\n";
		(*testout) << "connected = " << "\n" << connectedpairs << "\n";
		}



		// loch = CalcLocH (locpoints, locfaces, h);

		stat.nff = adfront.GetNF();
		stat.vol = adfront.Volume();
		if (stat.vol < 0)
		{
			break;
		}

		oldnp = locpoints.Size();
		oldnf = locfaces.Size();


		allowpoint.SetSize(locpoints.Size());
		if (uselocalh != 0 && stat.qualclass <= 3)
		{
	  for (int i = 1; i <= allowpoint.Size(); i++)
	  {
		  allowpoint.Elem(i) = (mesh.GetH(locpoints.Get(i)) > 0.4 * hshould / mp.sloppy) ? 2 : 1;
	  }
		}
		else
		{
	  allowpoint = 2;
		}



		if (stat.qualclass >= mp.starshapeclass && mp.baseelnp != 4)
		{
		NgProfiler.RegionTimer reg1 = new NgProfiler.RegionTimer(GenerateMesh_meshing3_timer_b);
		// star-shaped domain removing

		grouppoints.SetSize(0);
		groupfaces.SetSize(0);
		grouppindex.SetSize(0);
		groupfindex.SetSize(0);

		adfront.GetGroup(findex[0], grouppoints, groupfaces, grouppindex, groupfindex);

		bool onlytri = true;
		foreach (var i in groupfaces.Range())
		{
		  if (groupfaces[i].GetNP() != 3)
		  {
			onlytri = false;
		  }
		}

		if (onlytri && groupfaces.Size() <= 20 + 2 * stat.qualclass && netgen.GlobalMembers.FindInnerPoint(grouppoints, groupfaces, inp))
		{
			(*testout) << "inner point found" << "\n";

			for (int i = 1; i <= groupfaces.Size(); i++)
			{
		  adfront.DeleteFace(groupfindex.Get(i));
			}

			for (int i = 1; i <= groupfaces.Size(); i++)
			{
		  for (j = 1; j <= locfaces.Size(); j++)
		  {
			if (findex.Get(j) == groupfindex.Get(i))
			{
			  delfaces.Append(j);
			}
		  }
			}


			delfaces.SetSize(0);

			int npi;
			Element newel = new Element(TET);

			npi = mesh.AddPoint(inp);
			newel.SetNP(4);
			newel.PNum(4) = npi;

			for (int i = 1; i <= groupfaces.Size(); i++)
			{
			for (j = 1; j <= 3; j++)
			{
				newel.PNum(j) = adfront.GetGlobalIndex(grouppindex.Get(groupfaces.Get(i).PNum(j)));
			}
			mesh.AddVolumeElement(newel);
			}
			continue;
		}
		}

		found = 0;
		hasfound = 0;
		minerr = 1e6;

		//      int optother = 0;

		/*
		for(int i = 1; i <= locfaces.Size(); i++)
	  {
		(*testout) << "Face " << i << ": ";
		for (j = 1; j <= locfaces.Get(i).GetNP(); j++)
		  (*testout) << pindex.Get(locfaces.Get(i).PNum(j)) << " ";
		(*testout) << endl;
	  }
		for(int i = 1; i <= locpoints.Size(); i++)
	  {
		(*testout) << "p" << i
			   << ", gi = " << pindex.Get(i)
			   << " = " << locpoints.Get(i) << endl;
	  }
	  */

		GlobalMembers.minother = 1e10;
		GlobalMembers.minwithoutother = 1e10;

		bool impossible = true;

		for (rotind = 1; rotind <= locfaces[0].GetNP(); rotind++)
		{
		// set transformatino to reference coordinates

		if (locfaces[0].GetNP() == 3)
		{
			trans.Set(locpoints[locfaces[0].PNumMod(1 + rotind)], locpoints[locfaces[0].PNumMod(2 + rotind)], locpoints[locfaces[0].PNumMod(3 + rotind)], hshould);
		}
		else
		{
			trans.Set(locpoints[locfaces[0].PNumMod(1 + rotind)], locpoints[locfaces[0].PNumMod(2 + rotind)], locpoints[locfaces[0].PNumMod(4 + rotind)], hshould);
		}

		// trans.ToPlain (locpoints, plainpoints);

			plainpoints.SetSize(locpoints.Size());
			foreach (var i in locpoints.Range())
			{
			  trans.ToPlain(locpoints[i], plainpoints[i]);
			}

			foreach (var i in allowpoint.Range())
			{
			  if (plainpoints[i].Z() > 0)
			  {
				allowpoint[i] = false;
			  }
			}

		stat.cnttrials++;


		if (stat.cnttrials % 100 == 0)
		{
			(*testout) << "\n";
			for (int i = 1; i <= canuse.Size(); i++)
			{
		  (*testout) << foundmap.Get(i) << "/" << canuse.Get(i) << "/" << ruleused.Get(i) << " map/can/use rule " << rules.Get(i).Name() << "\n";
			}
			(*testout) << "\n";
		}

		// NgProfiler::StartTimer (meshing3_timer_c);
			GenerateMesh_meshing3_timer_c.Start();

		found = ApplyRules(plainpoints, allowpoint, locfaces, locfacesplit, connectedpairs, locelements, delfaces, stat.qualclass, mp.sloppy, rotind, err);

		if (found >= 0)
		{
			impossible = false;
		}
		if (found < 0)
		{
			found = 0;
		}

			GenerateMesh_meshing3_timer_c.Stop();
		// NgProfiler::StopTimer (meshing3_timer_c);

		if (found == 0)
		{
			loktestmode = false;
		}

		NgProfiler.RegionTimer reg2 = new NgProfiler.RegionTimer(GenerateMesh_meshing3_timer_d);

		if (loktestmode)
		{
			(*testout) << "plainpoints = " << "\n" << plainpoints << "\n";
			(*testout) << "Applyrules found " << found << "\n";
		}

		if (found != 0)
		{
			stat.cntsucc++;
		}

		locpoints.SetSize(plainpoints.Size());
		for (int i = oldnp + 1; i <= plainpoints.Size(); i++)
		{
		  trans.FromPlain(plainpoints.Elem(i), locpoints.Elem(i));
		}



		// avoid meshing from large to small mesh-size
		if (uselocalh != 0 && found != 0 && stat.qualclass <= 3)
		{
			for (int i = 1; i <= locelements.Size(); i++)
			{
			Point3d pmin = locpoints[locelements.Get(i).PNum(1)];
			Point3d pmax = new Point3d(pmin);
			for (j = 2; j <= 4; j++)
			{
				Point3d hp = locpoints[locelements.Get(i).PNum(j)];
				pmin.SetToMin(hp);
				pmax.SetToMax(hp);
			}

			if (mesh.GetMinH(pmin, pmax) < 0.4 * hshould / mp.sloppy)
			{
			  found = 0;
			}
			}
		}
		if (found != 0)
		{
			for (int i = 1; i <= locelements.Size(); i++)
			{
		  for (int j = 1; j <= 4; j++)
		  {
			  Point3d hp = locpoints[locelements.Get(i).PNum(j)];
			  if (netgen.GlobalMembers.Dist(hp, pmid) > hinner)
			  {
				found = 0;
			  }
		  }
			}
		}


		if (found != 0)
		{
		  ruleused.Elem(found)++;
		}


		// plotstat->Plot(stat);

		if (stat.cntelem != plotstat_oldne)
		{
			plotstat_oldne = stat.cntelem;

			PrintMessageCR(5, "El: ", stat.cntelem, " faces: ", stat.nff, " vol = ", (float)(100 * stat.vol / stat.vol0));

			multithread.percent = 100 - 100.0 * stat.vol / stat.vol0;
		}


		if (found != 0 && (hasfound == 0 || err < minerr))
		{

			if (testmode)
			{
			(*testout) << "found is active, 3" << "\n";
			for (int i = 1; i <= plainpoints.Size(); i++)
			{
				(*testout) << "p";
				if (i <= pindex.Size())
				{
			  (*testout) << pindex.Get(i) << ": ";
				}
				else
				{
			  (*testout) << "new: ";
				}
				(*testout) << plainpoints.Get(i) << "\n";
			}
			}



			hasfound = found;
			minerr = err;

			tempnewpoints.SetSize(0);
			for (int i = oldnp + 1; i <= locpoints.Size(); i++)
			{
		  tempnewpoints.Append(locpoints.Get(i));
			}

			tempnewfaces.SetSize(0);
			for (int i = oldnf + 1; i <= locfaces.Size(); i++)
			{
		  tempnewfaces.Append(locfaces.Get(i));
			}

			tempdelfaces.SetSize(0);
			for (int i = 1; i <= delfaces.Size(); i++)
			{
		  tempdelfaces.Append(delfaces.Get(i));
			}

			templocelements.SetSize(0);
			for (int i = 1; i <= locelements.Size(); i++)
			{
		  templocelements.Append(locelements.Get(i));
			}

			/*
			optother =
		  strcmp (problems[found], "other") == 0;
			*/
		}

		locpoints.SetSize(oldnp);
		locfaces.SetSize(oldnf);
		delfaces.SetSize(0);
		locelements.SetSize(0);
		}



		if (hasfound != 0)
		{

		/*
		if (optother)
		  (*testout) << "Other is optimal" << endl;
  
		if (minother < minwithoutother)
		  {
		    (*testout) << "Other is better, " << minother << " less " << minwithoutother << endl;
		  }
		  */

		for (int i = 1; i <= tempnewpoints.Size(); i++)
		{
		  locpoints.Append(tempnewpoints.Get(i));
		}
		for (int i = 1; i <= tempnewfaces.Size(); i++)
		{
		  locfaces.Append(tempnewfaces.Get(i));
		}
		for (int i = 1; i <= tempdelfaces.Size(); i++)
		{
		  delfaces.Append(tempdelfaces.Get(i));
		}
		for (int i = 1; i <= templocelements.Size(); i++)
		{
		  locelements.Append(templocelements.Get(i));
		}


		if (loktestmode)
		{
			(*testout) << "apply rule" << "\n";
			for (int i = 1; i <= locpoints.Size(); i++)
			{
			(*testout) << "p";
			if (i <= pindex.Size())
			{
			  (*testout) << pindex.Get(i) << ": ";
			}
			else
			{
			  (*testout) << "new: ";
			}
			(*testout) << locpoints.Get(i) << "\n";
			}
		}



		pindex.SetSize(locpoints.Size());

		for (int i = oldnp + 1; i <= locpoints.Size(); i++)
		{
			PointIndex globind = mesh.AddPoint(locpoints.Get(i));
			pindex.Elem(i) = adfront.AddPoint(locpoints.Get(i), globind);
		}

		for (int i = 1; i <= locelements.Size(); i++)
		{
			Point3d hp1;
			Point3d hp2;
			Point3d hp3;
			Point3d hp4;
			hp1 = locpoints[locelements.Get(i).PNum(1)];
			hp2 = locpoints[locelements.Get(i).PNum(2)];
			hp3 = locpoints[locelements.Get(i).PNum(3)];
			hp4 = locpoints[locelements.Get(i).PNum(4)];

			tetvol += (1.0 / 6.0) * (netgen.GlobalMembers.Cross(hp2 - hp1, hp3 - hp1) * (hp4 - hp1));

			for (j = 1; j <= locelements.Get(i).NP(); j++)
			{
		  locelements.Elem(i).PNum(j) = adfront.GetGlobalIndex(pindex[locelements.Get(i).PNum(j)]);
			}

			mesh.AddVolumeElement(locelements.Get(i));
			stat.cntelem++;
		}

		for (int i = oldnf + 1; i <= locfaces.Size(); i++)
		{
			for (j = 1; j <= locfaces.Get(i).GetNP(); j++)
			{
		  locfaces.Elem(i).PNum(j) = pindex[locfaces.Get(i).PNum(j)];
			}
			// (*testout) << "add face " << locfaces.Get(i) << endl;
			adfront.AddFace(locfaces.Get(i));
		}

		for (int i = 1; i <= delfaces.Size(); i++)
		{
		  adfront.DeleteFace(findex.Get(delfaces.Get(i)));
		}
		}
		else
		{
		adfront.IncrementClass(findex.Get(1));
		if (impossible && mp.check_impossible)
		{
			(*testout) << "skip face since it is impossible" << "\n";
			for (j = 0; j < 100; j++)
			{
		  adfront.IncrementClass(findex.Get(1));
			}
		}
		}

		locelements.SetSize(0);
		delpoints.SetSize(0);
		delfaces.SetSize(0);

		if (stat.qualclass >= mp.giveuptol)
		{
	  break;
		}
	}

	PrintMessage(5, ""); // line feed after statistics

	for (int i = 1; i <= ruleused.Size(); i++)
	{
	  (*testout) << setw(4) << ruleused.Get(i) << " times used rule " << rules.Get(i).Name() << "\n";
	}


	if (!mp.baseelnp && adfront.Empty())
	{
	  return MESHING3_OK;
	}

	if (mp.baseelnp && adfront.Empty(mp.baseelnp))
	{
	  return MESHING3_OK;
	}

	if (stat.vol < -1e-15)
	{
	  return MESHING3_NEGVOL;
	}

	return MESHING3_NEGVOL;
  }

  ///
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  int ApplyRules(Array<Point3d, PointIndex::BASE> lpoints, Array<int, PointIndex::BASE> allowpoint, Array<MiniElement2d> lfaces, INDEX lfacesplit, INDEX_2_HASHTABLE<int> connectedpairs, Array<Element> elements, Array<INDEX> delfaces, int tolerance, double sloppy, int rotind1, ref float retminerr);

  ///

  /*
    // was war das ????
  static double CalcLocH (const Array<Point3d> & locpoints,
  			const Array<MiniElement2d> & locfaces,
  			double h)
  {
    return h;
  
    
    int i, j;
    double hi, h1, d, dn, sum, weight, wi;
    Point3d p0, pc;
    Vec3d n, v1, v2;
  
    p0.X() = p0.Y() = p0.Z() = 0;
    for (j = 1; j <= 3; j++)
      {
        p0.X() += locpoints.Get(locfaces.Get(1).PNum(j)).X();
        p0.Y() += locpoints.Get(locfaces.Get(1).PNum(j)).Y();
        p0.Z() += locpoints.Get(locfaces.Get(1).PNum(j)).Z();
      }
    p0.X() /= 3; p0.Y() /= 3; p0.Z() /= 3;
    
    v1 = locpoints.Get(locfaces.Get(1).PNum(2)) -
      locpoints.Get(locfaces.Get(1).PNum(1));
    v2 = locpoints.Get(locfaces.Get(1).PNum(3)) -
      locpoints.Get(locfaces.Get(1).PNum(1));
  
    h1 = v1.Length();
    n = Cross (v2, v1);
    n /= n.Length();
  
    sum = 0;
    weight = 0;
  
    for(int i = 1; i <= locfaces.Size(); i++)
      {
        pc.X() = pc.Y() = pc.Z() = 0;
        for (j = 1; j <= 3; j++)
  	{
  	  pc.X() += locpoints.Get(locfaces.Get(i).PNum(j)).X();
  	  pc.Y() += locpoints.Get(locfaces.Get(i).PNum(j)).Y();
  	  pc.Z() += locpoints.Get(locfaces.Get(i).PNum(j)).Z();
  	}
        pc.X() /= 3; pc.Y() /= 3; pc.Z() /= 3;
  
        d = Dist (p0, pc);
        dn = n * (pc - p0);
        hi = Dist (locpoints.Get(locfaces.Get(i).PNum(1)),
  		 locpoints.Get(locfaces.Get(i).PNum(2)));
  		 
        if (dn > -0.2 * h1)
  	{
  	  wi = 1 / (h1 + d);
  	  wi *= wi;
  	}
        else
  	wi = 0;
  
        sum += hi * wi;
        weight += wi;
      }
  
    return sum/weight;
  }
  */

  public PointIndex AddPoint(Point3d p, PointIndex globind)
  {
	return adfront.AddPoint(p, globind);
  }

  ///
  public void AddBoundaryElement(Element2d elem)
  {
	MiniElement2d mini = new MiniElement2d(elem.GetNP());
	for (int j = 0; j < elem.GetNP(); j++)
	{
	  mini[j] = elem[j];
	}
	adfront.AddFace(mini);
  }

  ///
  public void AddBoundaryElement(MiniElement2d elem)
  {
	adfront.AddFace(elem);
  }

  ///
  public int AddConnectedPair(INDEX_2 apair)
  {
	return adfront.AddConnectedPair(apair);
  }

  ///
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
//C++ TO C# CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
  Timer BlockFill_t("Mesing3::BlockFill");
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private int[][] BlockFill_elind =
  {
	  new int[] {1, 8, 2, 4},
	  new int[] {1, 8, 4, 3},
	  new int[] {1, 8, 3, 7},
	  new int[] {1, 8, 7, 5},
	  new int[] {1, 8, 5, 6},
	  new int[] {1, 8, 6, 2}
  };

  public void BlockFill(Mesh mesh, double gh)
  {
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static Timer t("Mesing3::BlockFill");
	RegionTimer reg = new RegionTimer(BlockFill_t);

	PrintMessage(3, "Block-filling called (obsolete) ");

	int i;
	int j = 0;
	int i1;
	int i2;
	int i3;
	int j1;
	int j2;
	int j3;
	int n1;
	int n2;
	int n3;
	int n;
	int min1;
	int min2;
	int min3;
	int max1;
	int max2;
	int max3;
	int changed;
	int filled;
	double xmin = 0;
	double xmax = 0;
	double ymin = 0;
	double ymax = 0;
	double zmin = 0;
	double zmax = 0;
	double xminb;
	double xmaxb;
	double yminb;
	double ymaxb;
	double zminb;
	double zmaxb;
	//double rad = 0.7 * gh;

	for (int i = 1; i <= adfront.GetNP(); i++)
	{
		Point3d p = adfront.GetPoint(new PointIndex(i));
		if (i == 1)
		{
		xmin = xmax = p.X();
		ymin = ymax = p.Y();
		zmin = zmax = p.Z();
		}
		else
		{
		if (p.X() < xmin)
		{
			xmin = p.X();
		}
		if (p.X() > xmax)
		{
			xmax = p.X();
		}
		if (p.Y() < ymin)
		{
			ymin = p.Y();
		}
		if (p.Y() > ymax)
		{
			ymax = p.Y();
		}
		if (p.Z() < zmin)
		{
			zmin = p.Z();
		}
		if (p.Z() > zmax)
		{
			zmax = p.Z();
		}
		}
	}

	xmin -= 5 * gh;
	ymin -= 5 * gh;
	zmin -= 5 * gh;

	n1 = (int)((xmax - xmin) / gh + 5);
	n2 = (int)((ymax - ymin) / gh + 5);
	n3 = (int)((zmax - zmin) / gh + 5);
	n = n1 * n2 * n3;

	PrintMessage(5, "n1 = ", n1, " n2 = ", n2, " n3 = ", n3);

	Array<blocktyp> inner = new Array<blocktyp>(n);
	Array<PointIndex> pointnr = new Array<PointIndex>(n);
	Array<int> frontpointnr = new Array<int>(n);


	// initialize inner to 1

	for (int i = 1; i <= n; i++)
	{
	  inner.Elem(i) = BLOCKUNDEF;
	}


	// set blocks cutting surfaces to 0

	for (int i = 1; i <= adfront.GetNF(); i++)
	{
		MiniElement2d el = adfront.GetFace(i);
		xminb = xmax;
		xmaxb = xmin;
		yminb = ymax;
		ymaxb = ymin;
		zminb = zmax;
		zmaxb = zmin;

		for (j = 1; j <= 3; j++)
		{
		Point3d p = adfront.GetPoint(el.PNum(j));
		if (p.X() < xminb)
		{
			xminb = p.X();
		}
		if (p.X() > xmaxb)
		{
			xmaxb = p.X();
		}
		if (p.Y() < yminb)
		{
			yminb = p.Y();
		}
		if (p.Y() > ymaxb)
		{
			ymaxb = p.Y();
		}
		if (p.Z() < zminb)
		{
			zminb = p.Z();
		}
		if (p.Z() > zmaxb)
		{
			zmaxb = p.Z();
		}
		}



		double filldist = 0.2; // globflags.GetNumFlag ("filldist", 0.4);
		xminb -= filldist * gh;
		xmaxb += filldist * gh;
		yminb -= filldist * gh;
		ymaxb += filldist * gh;
		zminb -= filldist * gh;
		zmaxb += filldist * gh;

		min1 = (int)((xminb - xmin) / gh) + 1;
		max1 = (int)((xmaxb - xmin) / gh) + 1;
		netgen.GlobalMembers.min2 = (int)((yminb - ymin) / gh) + 1;
		netgen.GlobalMembers.max2 = (int)((ymaxb - ymin) / gh) + 1;
		netgen.GlobalMembers.min3 = (int)((zminb - zmin) / gh) + 1;
		netgen.GlobalMembers.max3 = (int)((zmaxb - zmin) / gh) + 1;


		for (i1 = min1; i1 <= max1; i1++)
		{
	  for (i2 = netgen.GlobalMembers.min2; i2 <= netgen.GlobalMembers.max2; i2++)
	  {
		for (i3 = netgen.GlobalMembers.min3; i3 <= netgen.GlobalMembers.max3; i3++)
		{
		  inner.Elem(i3 + (i2 - 1) * n3 + (i1 - 1) * n2 * n3) = BLOCKBOUND;
		}
	  }
		}
	}




	while (true)
	{
		int undefi = 0;
		Point3d undefp = new Point3d();

		for (i1 = 1; i1 <= n1 && !undefi; i1++)
		{
	  for (i2 = 1; i2 <= n2 && !undefi; i2++)
	  {
		for (i3 = 1; i3 <= n3 && !undefi; i3++)
		{
			i = i3 + (i2 - 1) * n3 + (i1 - 1) * n2 * n3;
			if (inner.Elem(i) == BLOCKUNDEF)
			{
			undefi = i;
			undefp.X() = xmin + (i1 - 0.5) * gh;
			undefp.Y() = ymin + (i2 - 0.5) * gh;
			undefp.Z() = zmin + (i3 - 0.5) * gh;
			}
		}
	  }
		}

		if (undefi == 0)
		{
	  break;
		}

		//      PrintMessage (5, "Test point: ", undefp);

		if (adfront.Inside(undefp))
		{
		//	  (*mycout) << "inner" << endl;
		inner.Elem(undefi) = BLOCKINNER;
		}
		else
		{
		//	  (*mycout) << "outer" << endl;
		inner.Elem(undefi) = BLOCKOUTER;
		}

		do
		{
		changed = 0;
		for (i1 = 1; i1 <= n1; i1++)
		{
		  for (i2 = 1; i2 <= n2; i2++)
		  {
			for (i3 = 1; i3 <= n3; i3++)
			{
			i = i3 + (i2 - 1) * n3 + (i1 - 1) * n2 * n3;

			for (int k = 1; k <= 3; k++)
			{
				switch (k)
				{
			  case 1:
				  j = i + n2 * n3;
				  break;
			  case 2:
				  j = i + n3;
				  break;
			  case 3:
				  j = i + 1;
				  break;
				}

				if (j > n1 * n2 * n3)
				{
					continue;
				}

				if (inner.Elem(i) == BLOCKOUTER && inner.Elem(j) == BLOCKUNDEF)
				{
				changed = 1;
				inner.Elem(j) = BLOCKOUTER;
				}
				if (inner.Elem(j) == BLOCKOUTER && inner.Elem(i) == BLOCKUNDEF)
				{
				changed = 1;
				inner.Elem(i) = BLOCKOUTER;
				}
				if (inner.Elem(i) == BLOCKINNER && inner.Elem(j) == BLOCKUNDEF)
				{
				changed = 1;
				inner.Elem(j) = BLOCKINNER;
				}
				if (inner.Elem(j) == BLOCKINNER && inner.Elem(i) == BLOCKUNDEF)
				{
				changed = 1;
				inner.Elem(i) = BLOCKINNER;
				}
			}
			}
		  }
		}
		} while (changed != 0);

	}



	filled = 0;
	for (int i = 1; i <= n; i++)
	{
	  if (inner.Elem(i) == BLOCKINNER)
	  {
	  filled++;
	  }
	}
	PrintMessage(5, "Filled blocks: ", filled);

	for (int i = 1; i <= n; i++)
	{
		pointnr.Elem(i) = 0;
		frontpointnr.Elem(i) = 0;
	}

	for (i1 = 1; i1 <= n1 - 1; i1++)
	{
	  for (i2 = 1; i2 <= n2 - 1; i2++)
	  {
		for (i3 = 1; i3 <= n3 - 1; i3++)
		{
		i = i3 + (i2 - 1) * n3 + (i1 - 1) * n2 * n3;
		if (inner.Elem(i) == BLOCKINNER)
		{
			for (j1 = i1; j1 <= i1 + 1; j1++)
			{
		  for (j2 = i2; j2 <= i2 + 1; j2++)
		  {
			for (j3 = i3; j3 <= i3 + 1; j3++)
			{
				j = j3 + (j2 - 1) * n3 + (j1 - 1) * n2 * n3;
				if (pointnr.Get(j) == 0)
				{
				Point3d hp = new Point3d(xmin + (j1 - 1) * gh, ymin + (j2 - 1) * gh, zmin + (j3 - 1) * gh);
				pointnr.Elem(j) = mesh.AddPoint(hp);
				frontpointnr.Elem(j) = AddPoint(hp, pointnr.Elem(j));

				}
			}
		  }
			}
		}
		}
	  }
	}


	for (i1 = 2; i1 <= n1 - 1; i1++)
	{
	  for (i2 = 2; i2 <= n2 - 1; i2++)
	  {
		for (i3 = 2; i3 <= n3 - 1; i3++)
		{
		i = i3 + (i2 - 1) * n3 + (i1 - 1) * n2 * n3;
		if (inner.Elem(i) == BLOCKINNER)
		{
			int[] pn = new int[9];
			pn[1] = pointnr.Get(i);
			pn[2] = pointnr.Get(i + 1);
			pn[3] = pointnr.Get(i + n3);
			pn[4] = pointnr.Get(i + n3 + 1);
			pn[5] = pointnr.Get(i + n2 * n3);
			pn[6] = pointnr.Get(i + n2 * n3 + 1);
			pn[7] = pointnr.Get(i + n2 * n3 + n3);
			pn[8] = pointnr.Get(i + n2 * n3 + n3 + 1);
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//			static int elind[][4] = { { 1, 8, 2, 4 }, { 1, 8, 4, 3 }, { 1, 8, 3, 7 }, { 1, 8, 7, 5 }, { 1, 8, 5, 6 }, { 1, 8, 6, 2 } };
			for (j = 1; j <= 6; j++)
			{
			Element el = new Element(4);
			for (int k = 1; k <= 4; k++)
			{
			  el.PNum(k) = pn[BlockFill_elind[j - 1][k - 1]];
			}

			mesh.AddVolumeElement(el);
			}
		}
		}
	  }
	}



	for (i1 = 2; i1 <= n1 - 1; i1++)
	{
	  for (i2 = 2; i2 <= n2 - 1; i2++)
	  {
		for (i3 = 2; i3 <= n3 - 1; i3++)
		{
		i = i3 + (i2 - 1) * n3 + (i1 - 1) * n2 * n3;
		if (inner.Elem(i) == BLOCKINNER)
		{
			int pi1 = 0;
			int pi2 = 0;
			int pi3 = 0;
			int pi4 = 0;

			int pn1 = frontpointnr.Get(i);
			int pn2 = frontpointnr.Get(i + 1);
			int pn3 = frontpointnr.Get(i + n3);
			int pn4 = frontpointnr.Get(i + n3 + 1);
			int pn5 = frontpointnr.Get(i + n2 * n3);
			int pn6 = frontpointnr.Get(i + n2 * n3 + 1);
			int pn7 = frontpointnr.Get(i + n2 * n3 + n3);
			int pn8 = frontpointnr.Get(i + n2 * n3 + n3 + 1);

			for (int k = 1; k <= 6; k++)
			{
			switch (k)
			{
			  case 1: // j3 = i3+1
				j = i + 1;
				pi1 = pn2;
				pi2 = pn6;
				pi3 = pn4;
				pi4 = pn8;
				break;
			  case 2: // j3 = i3-1
				j = i - 1;
				pi1 = pn1;
				pi2 = pn3;
				pi3 = pn5;
				pi4 = pn7;
				break;
			  case 3: // j2 = i2+1
				j = i + n3;
				pi1 = pn3;
				pi2 = pn4;
				pi3 = pn7;
				pi4 = pn8;
				break;
			  case 4: // j2 = i2-1
				j = i - n3;
				pi1 = pn1;
				pi2 = pn5;
				pi3 = pn2;
				pi4 = pn6;
				break;
			  case 5: // j1 = i1+1
				j = i + n3 * n2;
				pi1 = pn5;
				pi2 = pn7;
				pi3 = pn6;
				pi4 = pn8;
				break;
			  case 6: // j1 = i1-1
				j = i - n3 * n2;
				pi1 = pn1;
				pi2 = pn2;
				pi3 = pn3;
				pi4 = pn4;
				break;
			}

			if (inner.Get(j) == BLOCKBOUND)
			{
				MiniElement2d face = new MiniElement2d();
				face.PNum(1) = pi4;
				face.PNum(2) = pi1;
				face.PNum(3) = pi3;
				AddBoundaryElement(face);

				face.PNum(1) = pi1;
				face.PNum(2) = pi4;
				face.PNum(3) = pi2;
				AddBoundaryElement(face);

			}
			}
		}
		}
	  }
	}
  }

  ///

  /*
  static const AdFront3 * locadfront;
  static int TestInner (const Point3d & p)
  {
    return locadfront->Inside (p);
  }
  static int TestSameSide (const Point3d & p1, const Point3d & p2)
  {
    return locadfront->SameSide (p1, p2);
  }
  */



//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
//C++ TO C# CONVERTER TODO TASK: The following statement was not recognized, possibly due to an unrecognized macro:
  Timer BlockFillLocalH_t("Mesing3::BlockFillLocalH");

  public void BlockFillLocalH(Mesh mesh, MeshingParameters mp)
  {
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static Timer t("Mesing3::BlockFillLocalH");
	RegionTimer reg = new RegionTimer(BlockFillLocalH_t);

	double filldist = mp.filldist;

	// (*testout) << "blockfill local h" << endl;
	// (*testout) << "rel filldist = " << filldist << endl;
	PrintMessage(3, "blockfill local h");


	Array<Point < 3>> npoints = new Array<Point < 3>>();

	adfront.CreateTrees();

	Box < 3> bbox(Box < 3>.EMPTY_BOX);
	double maxh = 0;

	for (int i = 1; i <= adfront.GetNF(); i++)
	{
		MiniElement2d el = adfront.GetFace(i);
		for (int j = 1; j <= 3; j++)
		{
		Point3d p1 = adfront.GetPoint(el.PNumMod(j));
		Point3d p2 = adfront.GetPoint(el.PNumMod(j + 1));

		double hi = netgen.GlobalMembers.Dist(p1, p2);
		if (hi > maxh)
		{
			maxh = hi;
		}

		bbox.Add(p1);
		}
	}


	Point3d mpmin = bbox.PMin();
	Point3d mpmax = bbox.PMax();
	Point3d mpc = netgen.GlobalMembers.Center(mpmin, mpmax);
	double d = netgen.GlobalMembers.max3(mpmax.X() - mpmin.X(), mpmax.Y() - mpmin.Y(), mpmax.Z() - mpmin.Z()) / 2;
//C++ TO C# CONVERTER TODO TASK: The following line was determined to be a copy assignment (rather than a reference assignment) - this should be verified and a 'CopyFrom' method should be created:
//ORIGINAL LINE: mpmin = mpc - Vec3d(d, d, d);
	mpmin.CopyFrom(mpc - new Vec3d(d, d, d));
//C++ TO C# CONVERTER TODO TASK: The following line was determined to be a copy assignment (rather than a reference assignment) - this should be verified and a 'CopyFrom' method should be created:
//ORIGINAL LINE: mpmax = mpc + Vec3d(d, d, d);
	mpmax.CopyFrom(mpc + new Vec3d(d, d, d));
	Box3d meshbox = new Box3d(mpmin, mpmax);

	LocalH loch2 = new LocalH(mpmin, mpmax, 1);

	if (mp.maxh < maxh)
	{
		maxh = mp.maxh;
	}

	bool changed;
	do
	{
		mesh.LocalHFunction().ClearFlags();

		for (int i = 1; i <= adfront.GetNF(); i++)
		{
		MiniElement2d el = adfront.GetFace(i);

		Box < 3> bbox(adfront.GetPoint(el[0]));
		bbox.Add(adfront.GetPoint(el[1]));
		bbox.Add(adfront.GetPoint(el[2]));


		double filld = filldist * bbox.Diam();
		bbox.Increase(filld);

			  mesh.LocalHFunction().CutBoundary(bbox); // .PMin(), bbox.PMax());
		}

		//      locadfront = adfront;
		mesh.LocalHFunction().FindInnerBoxes(adfront, null);

		npoints.SetSize(0);
		mesh.LocalHFunction().GetInnerPoints(npoints);

		changed = false;
		for (int i = 1; i <= npoints.Size(); i++)
		{
		if (mesh.LocalHFunction().GetH(npoints.Get(i)) > 1.5 * maxh)
		{
			mesh.LocalHFunction().SetH(npoints.Get(i), maxh);
			changed = true;
		}
		}
	} while (changed);

	if (debugparam.slowchecks)
	{
	  (*testout) << "Blockfill with points: " << "\n";
	}
	for (int i = 1; i <= npoints.Size(); i++)
	{
		if (meshbox.IsIn(npoints.Get(i)))
		{
		PointIndex gpnum = mesh.AddPoint(npoints.Get(i));
		adfront.AddPoint(npoints.Get(i), gpnum);

		if (debugparam.slowchecks)
		{
			(*testout) << npoints.Get(i) << "\n";
			if (!adfront.Inside(npoints.Get(i)))
			{
			Console.Write("add outside point");
			Console.Write("\n");
			(*testout) << "outside" << "\n";
			}
		}

		}
	}



	// find outer points

	loch2.ClearFlags();

	for (int i = 1; i <= adfront.GetNF(); i++)
	{
		MiniElement2d el = adfront.GetFace(i);
		Point3d pmin = adfront.GetPoint(el.PNum(1));
		Point3d pmax = new Point3d(pmin);

		for (int j = 2; j <= 3; j++)
		{
		Point3d p = adfront.GetPoint(el.PNum(j));
		pmin.SetToMin(p);
		pmax.SetToMax(p);
		}

		loch2.SetH(netgen.GlobalMembers.Center(pmin, pmax), netgen.GlobalMembers.Dist(pmin, pmax));
	}

	for (int i = 1; i <= adfront.GetNF(); i++)
	{
		MiniElement2d el = adfront.GetFace(i);
		Point3d pmin = adfront.GetPoint(el.PNum(1));
		Point3d pmax = new Point3d(pmin);

		for (int j = 2; j <= 3; j++)
		{
		Point3d p = adfront.GetPoint(el.PNum(j));
		pmin.SetToMin(p);
		pmax.SetToMax(p);
		}

		double filld = filldist * netgen.GlobalMembers.Dist(pmin, pmax);
//C++ TO C# CONVERTER TODO TASK: The following line was determined to be a copy assignment (rather than a reference assignment) - this should be verified and a 'CopyFrom' method should be created:
//ORIGINAL LINE: pmin = pmin - Vec3d(filld, filld, filld);
		pmin.CopyFrom(pmin - new Vec3d(filld, filld, filld));
//C++ TO C# CONVERTER TODO TASK: The following line was determined to be a copy assignment (rather than a reference assignment) - this should be verified and a 'CopyFrom' method should be created:
//ORIGINAL LINE: pmax = pmax + Vec3d(filld, filld, filld);
		pmax.CopyFrom(pmax + new Vec3d(filld, filld, filld));
		// loch2.CutBoundary (pmin, pmax);
		loch2.CutBoundary(Box < 3> (pmin, pmax)); // pmin, pmax);
	}

	// locadfront = adfront;
	loch2.FindInnerBoxes(adfront, null);

	npoints.SetSize(0);
	loch2.GetOuterPoints(npoints);

	for (int i = 1; i <= npoints.Size(); i++)
	{
		if (meshbox.IsIn(npoints.Get(i)))
		{
		PointIndex gpnum = mesh.AddPoint(npoints.Get(i));
		adfront.AddPoint(npoints.Get(i), gpnum);
		}
	}
  }

  /// uses points of adfront, and puts new elements into mesh
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  void Delaunay(Mesh mesh, int domainnr, MeshingParameters mp);
  ///
//C++ TO C# CONVERTER TODO TASK: C# has no concept of a 'friend' class:
//  friend class PlotVolMesh;
  ///
//C++ TO C# CONVERTER TODO TASK: C# has no concept of a 'friend' function:
//ORIGINAL LINE: friend void TestRules();
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  void TestRules();
}




/// status of mesh generation
public class MeshingStat3d
{
  ///
  public MeshingStat3d()
  {
	cntsucc = cnttrials = cntelem = qualclass = 0;
	vol0 = h = 1;
	problemindex = 1;
  }

  ///
  public int cntsucc;
  ///
  public int cnttrials;
  ///
  public int cntelem;
  ///
  public int nff;
  ///
  public int qualclass;
  ///
  public double vol0;
  ///
  public double vol;
  ///
  public double h;
  ///
  public int problemindex;
}





/*
template <typename POINTArray, typename FACEArray>
extern int FindInnerPoint (POINTArray & grouppoints,
			   FACEArray & groupfaces,
			   Point3d & p);

*/

















namespace netgen
{




public enum blocktyp
{
	BLOCKUNDEF,
	BLOCKINNER,
	BLOCKBOUND,
	BLOCKOUTER
}

}
