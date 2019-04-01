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
/* File:   meshing2.hpp                                                   */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Okt. 95                                                    */
/**************************************************************************/



public enum MESHING2_RESULT
{
  MESHING2_OK = 0,
  MESHING2_GIVEUP = 1
}


/*
   
The basis class for 2D mesh generation. 
Has the method GenerateMesh

For surface mesh generation, or non-Euklidean meshing,
derive from Meshing2, and replace transformation.

*/

public class Meshing2
{
  /// the current advancing front
  private AdFront2 adfront;
  /// rules for mesh generation
  private Array<netrule> rules = new Array<netrule>();
  /// statistics
  private Array<int> ruleused = new Array<int>();
  private Array<int> canuse = new Array<int>();
  private Array<int> foundmap = new Array<int>();
  /// 
  private Box < 3> boundingbox;
  ///
  private double starttime;
  ///
  private double maxarea;

  private Vec3d ex = new Vec3d();
  private Vec3d ey = new Vec3d();
  private Point3d globp1 = new Point3d();

  ///

  // global variable for visualization
//   static Array<Point3d> locpoints;
//   static Array<int> legalpoints;
//   static Array<Point2d> plainpoints;
//   static Array<int> plainzones;
//   static Array<INDEX_2> loclines;
//   // static int geomtrig;
//   //static const char * rname;
//   static int cntelem, trials, nfaces;
//   static int oldnl;
//   static int qualclass;


  public Meshing2(MeshingParameters mp, Box < 3> aboundingbox)
  {
	boundingbox = aboundingbox;

	LoadRules(null, mp.quad);
	// LoadRules ("rules/quad.rls");
	// LoadRules ("rules/triangle.rls");

	adfront = new AdFront2(boundingbox);
	starttime = GetTime();

	maxarea = -1;
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
	  if (rules[i] != null)
	  {
		  rules[i].Dispose();
	  }
	}
  }

  /// Load rules, either from file, or compiled rules
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  void LoadRules(string filename, bool quad);

  /// 
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private int GenerateMesh_timer = NgProfiler.CreateTimer("surface meshing");
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private int GenerateMesh_timer1 = NgProfiler.CreateTimer("surface meshing1");
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private int GenerateMesh_timer2 = NgProfiler.CreateTimer("surface meshing2");
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private int GenerateMesh_timer3 = NgProfiler.CreateTimer("surface meshing3");
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private int GenerateMesh_ts1 = NgProfiler.CreateTimer("surface meshing start 1");
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private int GenerateMesh_ts2 = NgProfiler.CreateTimer("surface meshing start 2");
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private int GenerateMesh_ts3 = NgProfiler.CreateTimer("surface meshing start 3");
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private double GenerateMesh_maxviolate = 0;
//C++ TO C# CONVERTER NOTE: This was formerly a static local variable declaration (not allowed in C#):
  private int GenerateMesh_mtonnode = 0;

  public MESHING2_RESULT GenerateMesh(Mesh mesh, MeshingParameters mp, double gh, int facenr)
  {
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static int timer = NgProfiler::CreateTimer("surface meshing");

//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static int timer1 = NgProfiler::CreateTimer("surface meshing1");
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static int timer2 = NgProfiler::CreateTimer("surface meshing2");
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static int timer3 = NgProfiler::CreateTimer("surface meshing3");

//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static int ts1 = NgProfiler::CreateTimer("surface meshing start 1");
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static int ts2 = NgProfiler::CreateTimer("surface meshing start 2");
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//	static int ts3 = NgProfiler::CreateTimer("surface meshing start 3");


	NgProfiler.RegionTimer reg = new NgProfiler.RegionTimer(GenerateMesh_timer);

	NgProfiler.StartTimer(GenerateMesh_ts1);

	Array<int> pindex = new Array<int>();
	Array<int> lindex = new Array<int>();
	Array<int> delpoints = new Array<int>();
	Array<int> dellines = new Array<int>();

	Array<PointGeomInfo> upgeominfo = new Array<PointGeomInfo>(); // unique info
	Array<MultiPointGeomInfo> mpgeominfo = new Array<MultiPointGeomInfo>(); // multiple info

	Array<Element2d> locelements = new Array<Element2d>();

	int z1;
	int z2;
	int oldnp = -1;
	bool found;
	int rulenr = -1;
	Point < 3> p1, p2;

	PointGeomInfo blgeominfo1;
	PointGeomInfo blgeominfo2;

	bool morerisc;
	bool debugflag;

	double h;
	double his;
	double hshould;


	Array<Point3d> locpoints = new Array<Point3d>();
	Array<int> legalpoints = new Array<int>();
	Array<Point2d> plainpoints = new Array<Point2d>();
	Array<int> plainzones = new Array<int>();
	Array<INDEX_2> loclines = new Array<INDEX_2>();
	int cntelem = 0;
	int trials = 0;
	int nfaces = 0;
	int oldnl = 0;
	int qualclass;



	// test for 3d overlaps
	BoxTree < 3> surfeltree(boundingbox.PMin(), boundingbox.PMax());

	Array<int> intersecttrias = new Array<int>();
	Array<Point3d> critpoints = new Array<Point3d>();

	// test for doubled edges
	//INDEX_2_HASHTABLE<int> doubleedge(300000);


	testmode = 0;

	StartMesh();

	Array<Point2d> chartboundpoints = new Array<Point2d>();
	Array<Point3d> chartboundpoints3d = new Array<Point3d>();
	Array<INDEX_2> chartboundlines = new Array<INDEX_2>();

	// illegal points: points with more then 50 elements per node
	int maxlegalpoint = -1;
	int maxlegalline = -1;
	Array<int,PointIndex.BASE> trigsonnode = new Array<int,PointIndex.BASE>();
	Array<int,PointIndex.BASE> illegalpoint = new Array<int,PointIndex.BASE>();

	trigsonnode.SetSize(mesh.GetNP());
	illegalpoint.SetSize(mesh.GetNP());

	trigsonnode = 0;
	illegalpoint = 0;

	double totalarea = Area();
	double meshedarea = 0;


	// search tree for surface elements:
	/*
	for (sei = 0; sei < mesh.GetNSE(); sei++)
	  {
	const Element2d & sel = mesh[sei];

	if (sel.IsDeleted()) continue;

	if (sel.GetIndex() == facenr)
	  {
		Box<3> box;
		box.Set ( mesh[sel[0]] );
		box.Add ( mesh[sel[1]] );
		box.Add ( mesh[sel[2]] );
		surfeltree.Insert (box, sei);
	  }
	  }
	*/
	Array<SurfaceElementIndex> seia = new Array<SurfaceElementIndex>();
	mesh.GetSurfaceElementsOfFace(facenr, seia);
	for (int i = 0; i < seia.Size(); i++)
	{
	Element2d sel = mesh[seia[i]];

	if (sel.IsDeleted())
	{
		continue;
	}

	Box < 3> box;
	box.Set(mesh[sel[0]]);
	box.Add(mesh[sel[1]]);
	box.Add(mesh[sel[2]]);
	surfeltree.Insert(box, seia[i]);
	}

	NgProfiler.StopTimer(GenerateMesh_ts1);
	NgProfiler.StartTimer(GenerateMesh_ts2);

	if (totalarea > 0 || maxarea > 0)
	{
	  meshedarea = mesh.SurfaceArea();
	}
	  /*
	  for (SurfaceElementIndex sei = 0; sei < mesh.GetNSE(); sei++)
	{
	  const Element2d & sel = mesh[sei];
	  if (sel.IsDeleted()) continue;

	  double trigarea = Cross ( mesh[sel[1]]-mesh[sel[0]],
					mesh[sel[2]]-mesh[sel[0]] ).Length() / 2;


	  if (sel.GetNP() == 4)
		trigarea += Cross (Vec3d (mesh.Point (sel.PNum(1)),
					  mesh.Point (sel.PNum(3))),
				   Vec3d (mesh.Point (sel.PNum(1)),
					  mesh.Point (sel.PNum(4)))).Length() / 2;
	  meshedarea += trigarea;
	}
	  */
	  // cout << "meshedarea = " << meshedarea << " =?= "
	  // << mesh.SurfaceArea() << endl;

	NgProfiler.StopTimer(GenerateMesh_ts2);
	NgProfiler.StartTimer(GenerateMesh_ts3);

	string savetask = multithread.task;
	multithread.task = "Surface meshing";

	adfront.SetStartFront();


	int plotnexttrial = 999;

	double meshedarea_before = meshedarea;

	NgProfiler.StopTimer(GenerateMesh_ts3);

	while (!adfront.Empty() && !multithread.terminate)
	{
	NgProfiler.RegionTimer reg1 = new NgProfiler.RegionTimer(GenerateMesh_timer1);

	if (multithread.terminate)
	{
	  throw new Exception("Meshing stopped");
	}

	// known for STL meshing
	if (totalarea > 0)
	{
	  multithread.percent = 100 * meshedarea / totalarea;
	}
	/*
	  else
	  multithread.percent = 0;
	*/

	locpoints.SetSize(0);
	loclines.SetSize(0);
	pindex.SetSize(0);
	lindex.SetSize(0);
	delpoints.SetSize(0);
	dellines.SetSize(0);
	locelements.SetSize(0);



	// plot statistics
	if (trials > plotnexttrial)
	{
		PrintMessage(5, "faces = ", nfaces, " trials = ", trials, " elements = ", mesh.GetNSE(), " els/sec = ", (mesh.GetNSE() / (GetTime() - starttime + 0.0001)));
		plotnexttrial += 1000;
	}


	// unique-pgi, multi-pgi
	upgeominfo.SetSize(0);
	mpgeominfo.SetSize(0);


	nfaces = adfront.GetNFL();
	trials++;


	if (trials % 1000 == 0)
	{
		(*testout) << "\n";
		for (int i = 1; i <= canuse.Size(); i++)
		{
		(*testout) << foundmap.Get(i) << "/" << canuse.Get(i) << "/" << ruleused.Get(i) << " map/can/use rule " << rules.Get(i).Name() << "\n";
		}
		(*testout) << "\n";
	}


	int baselineindex = adfront.SelectBaseLine(p1, p2, blgeominfo1, blgeominfo2, qualclass);


	found = true;

	his = netgen.GlobalMembers.Dist(p1, p2);

	Point3d pmid = netgen.GlobalMembers.Center(p1, p2);
	hshould = CalcLocalH(pmid, mesh.GetH(pmid));
	if (gh < hshould)
	{
		hshould = gh;
	}

	mesh.RestrictLocalH(pmid, hshould);

	h = hshould;

	double hinner = (3 + qualclass) * netgen.GlobalMembers.max2(his, hshould);

	adfront.GetLocals(baselineindex, locpoints, mpgeominfo, loclines, pindex, lindex, 2 * hinner);


	NgProfiler.RegionTimer reg2 = new NgProfiler.RegionTimer(GenerateMesh_timer2);

	//(*testout) << "h for locals: " << 2*hinner << endl;


	//(*testout) << "locpoints " << locpoints << endl;

	if (qualclass > mp.giveuptol2d)
	{
		PrintMessage(3, "give up with qualclass ", qualclass);
		PrintMessage(3, "number of frontlines = ", adfront.GetNFL());
		// throw NgException ("Give up 2d meshing");
		break;
	}

	/*
	if (found && qualclass > 60)
	  {
	    found = 0;
	  }
	*/
	//      morerisc = ((qualclass > 20) && (qualclass % 2 == 1));
	//      morerisc = 1;
	morerisc = false;


	PointIndex gpi1 = adfront.GetGlobalIndex(pindex.Get(loclines[0].I1()));
	PointIndex gpi2 = adfront.GetGlobalIndex(pindex.Get(loclines[0].I2()));


	debugflag = (debugparam.haltsegment && (((debugparam.haltsegmentp1 == gpi1) && (debugparam.haltsegmentp2 == gpi2)) || ((debugparam.haltsegmentp1 == gpi2) && (debugparam.haltsegmentp2 == gpi1)))) || (debugparam.haltnode && ((debugparam.haltsegmentp1 == gpi1) || (debugparam.haltsegmentp2 == gpi1)));


	if (debugparam.haltface && debugparam.haltfacenr == facenr)
	{
		debugflag = true;
		Console.Write("set debugflag");
		Console.Write("\n");
	}

	if (debugparam.haltlargequalclass && qualclass > 50)
	{
	  debugflag = true;
	}

	// problem recognition !
	if (found && (gpi1 < illegalpoint.Size() + PointIndex.BASE) && (gpi2 < illegalpoint.Size() + PointIndex.BASE))
	{
		if (illegalpoint[gpi1] || illegalpoint[gpi2])
		{
		  found = false;
		}
	}


	Point2d p12d = new Point2d();
	Point2d p22d = new Point2d();

	if (found)
	{
		oldnp = locpoints.Size();
		oldnl = loclines.Size();

		if (debugflag)
		{
		  (*testout) << "define new transformation" << "\n";
		}

		DefineTransformation(p1, p2, blgeominfo1, blgeominfo2);

		plainpoints.SetSize(locpoints.Size());
		plainzones.SetSize(locpoints.Size());

		// (*testout) << endl;

		if (debugflag)
		{
		*testout << "3d->2d transformation" << "\n";
		*testout << "3d points: " << "\n" << locpoints << "\n";
		}

		for (int i = 1; i <= locpoints.Size(); i++)
		{
		// (*testout) << "pindex(i) = " << pindex[i-1] << endl;
		TransformToPlain(locpoints.Get(i), mpgeominfo.Get(i), plainpoints.Elem(i), h, plainzones.Elem(i));
		//		(*testout) << mpgeominfo.Get(i).GetPGI(1).u << " " << mpgeominfo.Get(i).GetPGI(1).v << " ";
		//		(*testout) << plainpoints.Get(i).X() << " " << plainpoints.Get(i).Y() << endl;
		//(*testout) << "transform " << locpoints.Get(i) << " to " << plainpoints.Get(i).X() << " " << plainpoints.Get(i).Y() << endl;
		}
		//	    (*testout) << endl << endl << endl;


		if (debugflag)
		{
		  *testout << "2d points: " << "\n" << plainpoints << "\n";
		}


		p12d = plainpoints.Get(1);
		p22d = plainpoints.Get(2);

		/*
		// last idea on friday
		plainzones.Elem(1) = 0;
		plainzones.Elem(2) = 0;
		*/


		/*
		// old netgen:
		for (i = 2; i <= loclines.Size(); i++)  // don't remove first line
		{
		z1 = plainzones.Get(loclines.Get(i).I1());
		z2 = plainzones.Get(loclines.Get(i).I2());
		
		if (z1 && z2 && (z1 != z2) || (z1 == -1) || (z2 == -1) )
		{
		loclines.DeleteElement(i);
		lindex.DeleteElement(i);
		oldnl--;
		i--;
		}
		}

		// 	  for (i = 1; i <= plainpoints.Size(); i++)
		// 	    if (plainzones.Elem(i) == -1)
		// 	      plainpoints.Elem(i) = Point2d (1e4, 1e4);
		*/



		for (int i = 2; i <= loclines.Size(); i++) // don't remove first line
		{
		// (*testout) << "loclines(i) = " << loclines.Get(i).I1() << " - " << loclines.Get(i).I2() << endl;
		z1 = plainzones.Get(loclines.Get(i).I1());
		z2 = plainzones.Get(loclines.Get(i).I2());


		// one inner point, one outer
		if ((z1 >= 0) != (z2 >= 0))
		{
			int innerp = (z1 >= 0) ? 1 : 2;
			if (IsLineVertexOnChart(locpoints.Get(loclines.Get(i).I1()), locpoints.Get(loclines.Get(i).I2()), innerp, adfront.GetLineGeomInfo(lindex.Get(i), innerp)))
			{
			  // pgeominfo.Get(loclines.Get(i).I(innerp))))

			if (!morerisc)
			{
				// use one end of line
				int pini;
				int pouti;
				Vec2d v = new Vec2d();

				pini = loclines.Get(i).I(innerp);
				pouti = loclines.Get(i).I(3 - innerp);

				Point2d pin = new Point2d(plainpoints.Get(pini));
				Point2d pout = new Point2d(plainpoints.Get(pouti));
				v = pout - pin;
				double len = v.Length();
				if (len <= 1e-6)
				{
				  (*testout) << "WARNING(js): inner-outer: short vector" << "\n";
				}
				else
				{
				  v /= len;
				}

				/*
				// don't elongate line towards base-line !!
				if (Vec2d (pin, p12d) * v > 0 &&
				Vec2d (pin, p22d) * v > 0)
				v *= -1;
				*/

				Point2d newpout = pin + 1000 * v;
//C++ TO C# CONVERTER TODO TASK: The following line was determined to be a copy assignment (rather than a reference assignment) - this should be verified and a 'CopyFrom' method should be created:
//ORIGINAL LINE: newpout = pout;
				newpout.CopyFrom(pout);


				plainpoints.Append(newpout);
				Point3d pout3d = locpoints.Get(pouti);
				locpoints.Append(pout3d);

				plainzones.Append(0);
				pindex.Append(-1);
				oldnp++;
				loclines.Elem(i).I(3 - innerp) = oldnp;
			}
			else
			{
			  plainzones.Elem(loclines.Get(i).I(3 - innerp)) = 0;
			}


			//		  (*testout) << "inner - outer correction" << endl;
			}
			else
			{
			// remove line
			loclines.DeleteElement(i);
			lindex.DeleteElement(i);
			oldnl--;
			i--;
			}
		}

		else if ((z1 > 0 && z2 > 0 && (z1 != z2)) || ((z1 < 0) && (z2 < 0)))
		{
			loclines.DeleteElement(i);
			lindex.DeleteElement(i);
			oldnl--;
			i--;
		}
		}





		legalpoints.SetSize(plainpoints.Size());
		for (int i = 1; i <= legalpoints.Size(); i++)
		{
		  legalpoints.Elem(i) = 1;
		}

		double avy = 0;
		for (int i = 1; i <= plainpoints.Size(); i++)
		{
		  avy += plainpoints.Elem(i).Y();
		}
		avy *= 1.0 / plainpoints.Size();


		for (int i = 1; i <= plainpoints.Size(); i++)
		{
		if (plainzones.Elem(i) < 0)
		{
			plainpoints.Elem(i) = new Point2d(1e4, 1e4);
			legalpoints.Elem(i) = 0;
		}
		if (pindex.Elem(i) == -1)
		{
			legalpoints.Elem(i) = 0;
		}


		if (plainpoints.Elem(i).Y() < -1e-10 * avy) // changed
		{
			legalpoints.Elem(i) = 0;
		}
		}
		/*
		  for (i = 3; i <= plainpoints.Size(); i++)
		  if (sqr (plainpoints.Get(i).X()) + sqr (plainpoints.Get(i).Y())
		  > sqr (2 + 0.2 * qualclass))
		  legalpoints.Elem(i) = 0;
		*/

		/*
		 int clp = 0;
		 for (i = 1; i <= plainpoints.Size(); i++)
		 if (legalpoints.Get(i))
		 clp++;
		 (*testout) << "legalpts: " << clp << "/" << plainpoints.Size() << endl;

		 // sort legal/illegal lines
		 int lastleg = 2;
		 int firstilleg = oldnl;

		 while (lastleg < firstilleg)
		 {
		 while (legalpoints.Get(loclines.Get(lastleg).I1()) &&
		 legalpoints.Get(loclines.Get(lastleg).I2()) &&
		 lastleg < firstilleg)
		 lastleg++;
		 while ( ( !legalpoints.Get(loclines.Get(firstilleg).I1()) ||
		 !legalpoints.Get(loclines.Get(firstilleg).I2())) &&
		 lastleg < firstilleg)
		 firstilleg--;
		
		 if (lastleg < firstilleg)
		 {
		 swap (loclines.Elem(lastleg), loclines.Elem(firstilleg));
		 swap (lindex.Elem(lastleg), lindex.Elem(firstilleg));
		 }
		 }

		 (*testout) << "leglines " << lastleg << "/" << oldnl << endl;
		*/


		GetChartBoundary(chartboundpoints, chartboundpoints3d, chartboundlines, h);

		oldnp = plainpoints.Size();

		maxlegalpoint = locpoints.Size();
		maxlegalline = loclines.Size();



		if (mp.checkchartboundary)
		{
		for (int i = 1; i <= chartboundpoints.Size(); i++)
		{
			plainpoints.Append(chartboundpoints.Get(i));
			locpoints.Append(chartboundpoints3d.Get(i));
			legalpoints.Append(0);
		}


		for (int i = 1; i <= chartboundlines.Size(); i++)
		{
			INDEX_2 line = new INDEX_2(chartboundlines.Get(i).I1() + oldnp, chartboundlines.Get(i).I2() + oldnp);
			loclines.Append(line);
			//	      (*testout) << "line: " << line.I1() << "-" << line.I2() << endl;
		}
		}

		oldnl = loclines.Size();
		oldnp = plainpoints.Size();
	}


	/*
	  if (qualclass > 100)
	  {
	  multithread.drawing = 1;
	  glrender(1);
	  cout << "qualclass 100, nfl = " << adfront->GetNFL() << endl;
	  }
	*/

	if (found)
	{
		rulenr = ApplyRules(plainpoints, legalpoints, maxlegalpoint, loclines, maxlegalline, locelements, dellines, qualclass, mp);

		//	    (*testout) << "Rule Nr = " << rulenr << endl;
		if (rulenr == 0)
		{
		found = false;
		if (debugflag || debugparam.haltnosuccess)
		{
		  PrintWarning("no rule found");
		}
		}
	}

	NgProfiler.RegionTimer reg3 = new NgProfiler.RegionTimer(GenerateMesh_timer3);


	for (int i = 1; i <= locelements.Size() && found; i++)
	{
		Element2d el = locelements.Get(i);

		for (int j = 1; j <= el.GetNP(); j++)
		{
		  if (el.PNum(j) <= oldnp && pindex.Get(el.PNum(j)) == -1)
		  {
		  found = false;
		  PrintSysError("meshing2, index missing");
		  }
		}
	}


	if (found)
	{
		locpoints.SetSize(plainpoints.Size());
		upgeominfo.SetSize(locpoints.Size());

		for (int i = oldnp + 1; i <= plainpoints.Size(); i++)
		{
		int err = TransformFromPlain(plainpoints.Elem(i), locpoints.Elem(i), upgeominfo.Elem(i), h);

		if (err != 0)
		{
			found = false;

			if (debugflag || debugparam.haltnosuccess)
			{
			  PrintSysError("meshing2, Backtransformation failed");
			}

			break;
		}
		}
	}


	//      for (i = 1; i <= oldnl; i++)
	//        adfront -> ResetClass (lindex[i]);


	/*
	  double violateminh;
	  if (qualclass <= 10)
	  violateminh = 3;
	  else
	  violateminh = 3 * qualclass;

	  if (uselocalh && found) //  && qualclass <= 10)
	  {
	  for (i = 1; i <= locelements.Size(); i++)
	  {
	  Point3d pmin = locpoints.Get(locelements.Get(i).PNum(1));
	  Point3d pmax = pmin;
	  for (j = 2; j <= 3; j++)
	  {
	  const Point3d & hp =
	  locpoints.Get(locelements.Get(i).PNum(j));
	  pmin.SetToMin (hp);
	  pmax.SetToMax (hp);
	  }
	  double minh = mesh.GetMinH (pmin, pmax);
	  if (h > violateminh * minh)
	  {
	  found = 0;
	  loclines.SetSize (oldnl);
	  locpoints.SetSize (oldnp);
	  }
	  }
	  }
	*/


	if (found)
	{
		double violateminh = 3 + 0.1 * netgen.GlobalMembers.sqr(qualclass);
		double minh = 1e8;
		double newedgemaxh = 0;
		for (int i = oldnl + 1; i <= loclines.Size(); i++)
		{
		double eh = netgen.GlobalMembers.Dist(locpoints.Get(loclines.Get(i).I1()), locpoints.Get(loclines.Get(i).I2()));

		// Markus (brute force method to avoid bad elements on geometries like \_/ )
		//if(eh > 4.*mesh.GetH(locpoints.Get(loclines.Get(i).I1()))) found = 0;
		//if(eh > 4.*mesh.GetH(locpoints.Get(loclines.Get(i).I2()))) found = 0;
		// Markus end

		if (eh > newedgemaxh)
		{
		  newedgemaxh = eh;
		}
		}

		for (int i = 1; i <= locelements.Size(); i++)
		{
		Point3d pmin = locpoints.Get(locelements.Get(i).PNum(1));
		Point3d pmax = new Point3d(pmin);
		for (int j = 2; j <= locelements.Get(i).GetNP(); j++)
		{
			Point3d hp = locpoints.Get(locelements.Get(i).PNum(j));
			pmin.SetToMin(hp);
			pmax.SetToMax(hp);
		}
		double eh = mesh.GetMinH(pmin, pmax);
		if (eh < minh)
		{
		  minh = eh;
		}
		}

		for (int i = 1; i <= locelements.Size(); i++)
		{
		  for (int j = 1; j <= locelements.Get(i).GetNP(); j++)
		  {
		if (netgen.GlobalMembers.Dist2(locpoints.Get(locelements.Get(i).PNum(j)), pmid) > hinner * hinner)
		{
		  found = false;
		}
		  }
		}

		//	  cout << "violate = " << newedgemaxh / minh << endl;
//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//		static double maxviolate = 0;
		if (newedgemaxh / minh > GenerateMesh_maxviolate)
		{
		GenerateMesh_maxviolate = newedgemaxh / minh;
		//	      cout << "max minhviolate = " << maxviolate << endl;
		}


		if (newedgemaxh > violateminh * minh)
		{
		found = false;
		loclines.SetSize(oldnl);
		locpoints.SetSize(oldnp);

		if (debugflag || debugparam.haltnosuccess)
		{
		  PrintSysError("meshing2, maxh too large");
		}


		}
	}



	/*
	// test good ComputeLineGeoInfo
	if (found)
	{
	// is line on chart ?
	for (i = oldnl+1; i <= loclines.Size(); i++)
	{
	int gisize;
	void *geominfo;

	if (ComputeLineGeoInfo (locpoints.Get(loclines.Get(i).I1()),
	locpoints.Get(loclines.Get(i).I2()),
	gisize, geominfo))
	found = 0;
	}
	}
	*/


	// changed for OCC meshing
	if (found)
	{
		// take geominfo from dellines
		// upgeominfo.SetSize(locpoints.Size());

		/*
		  for (i = 1; i <= dellines.Size(); i++)
		  for (j = 1; j <= 2; j++)
		  {
		  upgeominfo.Elem(loclines.Get(dellines.Get(i)).I(j)) =
		  adfront -> GetLineGeomInfo (lindex.Get(dellines.Get(i)), j);
		  }
		*/


		for (int i = 1; i <= locelements.Size(); i++)
		{
		  for (int j = 1; j <= locelements.Get(i).GetNP(); j++)
		  {
		  int pi = locelements.Get(i).PNum(j);
		  if (pi <= oldnp)
		  {

			  if (ChooseChartPointGeomInfo(mpgeominfo.Get(pi), upgeominfo.Elem(pi)))
			  {
			  // cannot select, compute new one
			  PrintWarning("calc point geominfo instead of using");
			  if (ComputePointGeomInfo(locpoints.Get(pi), upgeominfo.Elem(pi)))
			  {
				  found = false;
				  PrintSysError("meshing2d, geominfo failed");
			  }
			  }
		  }
		  }
		}

		/*
		// use upgeominfo from ProjectFromPlane
		for (i = oldnp+1; i <= locpoints.Size(); i++)
		{
		if (ComputePointGeomInfo (locpoints.Get(i), upgeominfo.Elem(i)))
		{
		found = 0;
		if ( debugflag || debugparam.haltnosuccess )
		PrintSysError ("meshing2d, compute geominfo failed");
		}
		}
		*/
	}


	if (found && mp.checkoverlap)
	{
		// cout << "checkoverlap" << endl;
		// test for overlaps

		Point3d hullmin = new Point3d(1e10, 1e10, 1e10);
		Point3d hullmax = new Point3d(-1e10, -1e10, -1e10);

		for (int i = 1; i <= locelements.Size(); i++)
		{
		  for (int j = 1; j <= locelements.Get(i).GetNP(); j++)
		  {
		  Point3d p = locpoints.Get(locelements.Get(i).PNum(j));
		  hullmin.SetToMin(p);
		  hullmax.SetToMax(p);
		  }
		}
		hullmin += new Vec3d(-his, -his, -his);
		hullmax += new Vec3d(his, his, his);

		surfeltree.GetIntersecting(hullmin, hullmax, intersecttrias);

		critpoints.SetSize(0);
		for (int i = oldnp + 1; i <= locpoints.Size(); i++)
		{
		  critpoints.Append(locpoints.Get(i));
		}

		for (int i = 1; i <= locelements.Size(); i++)
		{
		Element2d tri = locelements.Get(i);
		if (tri.GetNP() == 3)
		{
			Point3d tp1 = locpoints.Get(tri.PNum(1));
			Point3d tp2 = locpoints.Get(tri.PNum(2));
			Point3d tp3 = locpoints.Get(tri.PNum(3));

			Vec3d tv1 = new Vec3d(tp1, tp2);
			Vec3d tv2 = new Vec3d(tp1, tp3);

			double lam1;
			double lam2;
			for (lam1 = 0.2; lam1 <= 0.8; lam1 += 0.2)
			{
			  for (lam2 = 0.2; lam2 + lam1 <= 0.8; lam2 += 0.2)
			  {
			  Point3d hp = tp1 + lam1 * tv1 + lam2 * tv2;
			  critpoints.Append(hp);
			  }
			}
		}
		else if (tri.GetNP() == 4)
		{
			Point3d tp1 = locpoints.Get(tri.PNum(1));
			Point3d tp2 = locpoints.Get(tri.PNum(2));
			Point3d tp3 = locpoints.Get(tri.PNum(3));
			Point3d tp4 = locpoints.Get(tri.PNum(4));

			double l1;
			double l2;
			for (l1 = 0.1; l1 <= 0.9; l1 += 0.1)
			{
			  for (l2 = 0.1; l2 <= 0.9; l2 += 0.1)
			  {
			  Point3d hp = new Point3d();
			  hp.X() = (1 - l1) * (1 - l2) * tp1.X() + l1 * (1 - l2) * tp2.X() + l1 * l2 * tp3.X() + (1 - l1) * l2 * tp4.X();
			  hp.Y() = (1 - l1) * (1 - l2) * tp1.Y() + l1 * (1 - l2) * tp2.Y() + l1 * l2 * tp3.Y() + (1 - l1) * l2 * tp4.Y();
			  hp.Z() = (1 - l1) * (1 - l2) * tp1.Z() + l1 * (1 - l2) * tp2.Z() + l1 * l2 * tp3.Z() + (1 - l1) * l2 * tp4.Z();


			  critpoints.Append(hp);
			  }
			}
		}
		}
		/*
		  for (i = oldnl+1; i <= loclines.Size(); i++)
		  {
		  Point3d hp = locpoints.Get(loclines.Get(i).I1());
		  Vec3d hv(hp, locpoints.Get(loclines.Get(i).I2()));
		  int ncp = 2;
		  for (j = 1; j <= ncp; j++)
		  critpoints.Append ( hp + (double(j)/(ncp+1)) * hv);
		  }
		*/


		/*
		  for (i = oldnp+1; i <= locpoints.Size(); i++)
		  {
		  const Point3d & p = locpoints.Get(i);
		*/


		for (int i = 1; i <= critpoints.Size(); i++)
		{
		Point3d p = critpoints.Get(i);

		for (int jj = 0; jj < intersecttrias.Size(); jj++)
		{
			// int j = intersecttrias.Get(jj);
			// const Element2d & el = mesh.SurfaceElement(j);

			SurfaceElementIndex j = intersecttrias[jj];
			Element2d el = mesh[j];

			int ntrig = (el.GetNP() == 3) ? 1 : 2;

			int jl;
			for (jl = 1; jl <= ntrig; jl++)
			{
			Point3d tp1 = new Point3d();
			Point3d tp2 = new Point3d();
			Point3d tp3 = new Point3d();

			if (jl == 1)
			{
				tp1 = new mesh.Point(el.PNum(1));
				tp2 = new mesh.Point(el.PNum(2));
				tp3 = new mesh.Point(el.PNum(3));
			}
			else
			{
				tp1 = new mesh.Point(el.PNum(1));
				tp2 = new mesh.Point(el.PNum(3));
				tp3 = new mesh.Point(el.PNum(4));
			}

			int onchart = 0;
			for (int k = 1; k <= el.GetNP(); k++)
			{
			  if (BelongsToActiveChart(new mesh.Point(el.PNum(k)), el.GeomInfoPi(k)))
			  {
				onchart = 1;
			  }
			}
			if (onchart == 0)
			{
			  continue;
			}

			Vec3d e1 = new Vec3d(tp1, tp2);
			Vec3d e2 = new Vec3d(tp1, tp3);
			Vec3d n = netgen.GlobalMembers.Cross(e1, e2);
			n /= n.Length();
			double lam1;
			double lam2;
			double lam3;
			lam3 = n * new Vec3d(tp1, p);
			LocalCoordinates(e1, e2, new Vec3d(tp1, p), lam1, lam2);

			if (ngsimd.GlobalMembers.fabs(lam3) < 0.1 * hshould && lam1 > 0 && lam2 > 0 && (lam1 + lam2) < 1)
			{
#if DEVELOP
				Console.Write("overlap");
				Console.Write("\n");
				(*testout) << "overlap:" << "\n" << "tri = " << tp1 << "-" << tp2 << "-" << tp3 << "\n" << "point = " << p << "\n" << "lam1, 2 = " << lam1 << ", " << lam2 << "\n" << "lam3 = " << lam3 << "\n";

				//		      cout << "overlap !!!" << endl;
#endif
				for (int k = 1; k <= 5; k++)
				{
				  adfront.IncrementClass(lindex.Get(1));
				}

				found = false;

				if (debugflag || debugparam.haltnosuccess)
				{
				  PrintWarning("overlapping");
				}


				if (debugparam.haltoverlap)
				{
				debugflag = true;
				}

				/*
				  multithread.drawing = 1;
				  glrender(1);
				*/
			}
			}
		}
		}
	}


	if (found)
	{
		// check, whether new front line already exists

		for (int i = oldnl + 1; i <= loclines.Size(); i++)
		{
		int nlgpi1 = loclines.Get(i).I1();
		int nlgpi2 = loclines.Get(i).I2();
		if (nlgpi1 <= pindex.Size() && nlgpi2 <= pindex.Size())
		{
			nlgpi1 = adfront.GetGlobalIndex(pindex.Get(nlgpi1));
			nlgpi2 = adfront.GetGlobalIndex(pindex.Get(nlgpi2));

			int exval = adfront.ExistsLine(nlgpi1, nlgpi2);
			if (exval != 0)
			{
			Console.Write("ERROR: new line exits, val = ");
			Console.Write(exval);
			Console.Write("\n");
			(*testout) << "ERROR: new line exits, val = " << exval << "\n";
			found = false;


			if (debugparam.haltexistingline)
			{
			  debugflag = true;
			}

			}
		}
		}

	}


	/*
	  if (found)
	  {
	  // check, whether new triangles insert edges twice
	  for (i = 1; i <= locelements.Size(); i++)
	  for (j = 1; j <= 3; j++)
	  {
	  int tpi1 = locelements.Get(i).PNumMod (j);
	  int tpi2 = locelements.Get(i).PNumMod (j+1);
	  if (tpi1 <= pindex.Size() && tpi2 <= pindex.Size())
	  {
	  tpi1 = adfront->GetGlobalIndex (pindex.Get(tpi1));
	  tpi2 = adfront->GetGlobalIndex (pindex.Get(tpi2));

	  if (doubleedge.Used (INDEX_2(tpi1, tpi2)))
	  {
	  if (debugparam.haltexistingline)
	  debugflag = 1;
	  cerr << "ERROR Insert edge "
	  << tpi1 << " - " << tpi2 << " twice !!!" << endl;
	  found = 0;
	  }
	  doubleedge.Set (INDEX_2(tpi1, tpi2), 1);
	  }
	  }
	  }
	*/


	if (found)
	{
		// everything is ok, perform mesh update

		ruleused.Elem(rulenr)++;


		pindex.SetSize(locpoints.Size());

		for (int i = oldnp + 1; i <= locpoints.Size(); i++)
		{
		PointIndex globind = mesh.AddPoint(locpoints.Get(i));
		pindex.Elem(i) = adfront.AddPoint(locpoints.Get(i), globind);
		}

		for (int i = oldnl + 1; i <= loclines.Size(); i++)
		{
		/*
		  for (j = 1; j <= locpoints.Size(); j++)
		  {
		  (*testout) << j << ": " << locpoints.Get(j) << endl;
		  }
		*/

		/*
		  ComputeLineGeoInfo (locpoints.Get(loclines.Get(i).I1()),
		  locpoints.Get(loclines.Get(i).I2()),
		  gisize, geominfo);
		*/

		if (pindex.Get(loclines.Get(i).I1()) == -1 || pindex.Get(loclines.Get(i).I2()) == -1)
		{
			(*testout) << "pindex is 0" << "\n";
		}

		if (!upgeominfo.Get(loclines.Get(i).I1()).trignum || !upgeominfo.Get(loclines.Get(i).I2()).trignum)
		{
			Console.Write("new el: illegal geominfo");
			Console.Write("\n");
		}

		adfront.AddLine(pindex.Get(loclines.Get(i).I1()), pindex.Get(loclines.Get(i).I2()), upgeominfo.Get(loclines.Get(i).I1()), upgeominfo.Get(loclines.Get(i).I2()));
		}
		for (int i = 1; i <= locelements.Size(); i++)
		{
		Element2d mtri = new Element2d(locelements.Get(i).GetNP());
		mtri = locelements.Get(i);
		mtri.SetIndex(facenr);


		// compute triangle geominfo:
		//	      (*testout) << "triggeominfo: ";
		for (int j = 1; j <= locelements.Get(i).GetNP(); j++)
		{
			mtri.GeomInfoPi(j) = upgeominfo.Get(locelements.Get(i).PNum(j));
			//		  (*testout) << mtri.GeomInfoPi(j).trignum << " ";
		}
		//	      (*testout) << endl;

		for (int j = 1; j <= locelements.Get(i).GetNP(); j++)
		{
			mtri.PNum(j) = locelements.Elem(i).PNum(j) = adfront.GetGlobalIndex(pindex.Get(locelements.Get(i).PNum(j)));
		}




		mesh.AddSurfaceElement(mtri);
		cntelem++;
		//	      cout << "elements: " << cntelem << endl;



		Box < 3> box;
		box.Set(mesh[mtri[0]]);
		box.Add(mesh[mtri[1]]);
		box.Add(mesh[mtri[2]]);
		surfeltree.Insert(box, mesh.GetNSE() - 1);

		Point3d sep1 = new mesh.Point(mtri.PNum(1));
		Point3d sep2 = new mesh.Point(mtri.PNum(2));
		Point3d sep3 = new mesh.Point(mtri.PNum(3));

		double trigarea = netgen.GlobalMembers.Cross(new Vec3d(sep1, sep2), new Vec3d(sep1, sep3)).Length() / 2;

		if (mtri.GetNP() == 4)
		{
			Point3d sep4 = new mesh.Point(mtri.PNum(4));
			trigarea += netgen.GlobalMembers.Cross(new Vec3d(sep1, sep3), new Vec3d(sep1, sep4)).Length() / 2;
		}

		meshedarea += trigarea;

		if (maxarea > 0 && meshedarea - meshedarea_before > maxarea)
		{
			cerr << "meshed area = " << meshedarea - meshedarea_before << "\n" << "maximal area = " << maxarea << "\n" << "GIVING UP" << "\n";
			return MESHING2_GIVEUP;
		}



		for (int j = 1; j <= locelements.Get(i).GetNP(); j++)
		{
			int gpi = locelements.Get(i).PNum(j);

			int oldts = trigsonnode.Size();
			if (gpi >= oldts + PointIndex.BASE)
			{
			trigsonnode.SetSize(gpi + 1 - PointIndex.BASE);
			illegalpoint.SetSize(gpi + 1 - PointIndex.BASE);
			for (int k = oldts + PointIndex.BASE; k <= gpi; k++)
			{
				trigsonnode[k] = 0;
				illegalpoint[k] = 0;
			}
			}

			trigsonnode[gpi]++;

			if (trigsonnode[gpi] > 20)
			{
			illegalpoint[gpi] = 1;
			//		      cout << "illegal point: " << gpi << endl;
			(*testout) << "illegal point: " << gpi << "\n";
			}

//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
//			static int mtonnode = 0;
			if (trigsonnode[gpi] > GenerateMesh_mtonnode)
			{
			  GenerateMesh_mtonnode = trigsonnode[gpi];
			}
		}
		//	      cout << "els = " << cntelem << " trials = " << trials << endl;
		//	      if (trials > 100)		return;
		}

		for (int i = 1; i <= dellines.Size(); i++)
		{
		  adfront.DeleteLine(lindex.Get(dellines.Get(i)));
		}

		//	  rname = rules.Get(rulenr)->Name();
#if MYGRAPH
		if (silentflag < 3)
		{
		plotsurf.DrawPnL(locpoints, loclines);
		plotsurf.Plot(testmode, testmode);
		}
#endif

		if (morerisc)
		{
		Console.Write("generated due to morerisc");
		Console.Write("\n");
		//	      multithread.drawing = 1;
		//	      glrender(1);
		}




		if (debugparam.haltsuccess || debugflag)
		{
		// adfront -> PrintOpenSegments (*testout);
		Console.Write("success of rule");
		Console.Write(rules.Get(rulenr).Name());
		Console.Write("\n");
		multithread.drawing = 1;
		multithread.testmode = 1;
		multithread.pause = 1;


		/*
		  extern STLGeometry * stlgeometry;
		  stlgeometry->ClearMarkedSegs();
		  for (i = 1; i <= loclines.Size(); i++)
		  {
		  stlgeometry->AddMarkedSeg(locpoints.Get(loclines.Get(i).I1()),
		  locpoints.Get(loclines.Get(i).I2()));
		  }
		*/

		(*testout) << "success of rule" << rules.Get(rulenr).Name() << "\n";
		(*testout) << "trials = " << trials << "\n";

		(*testout) << "locpoints " << "\n";
		for (int i = 1; i <= pindex.Size(); i++)
		{
		  (*testout) << adfront.GetGlobalIndex(pindex.Get(i)) << "\n";
		}

		(*testout) << "old number of lines = " << oldnl << "\n";
		for (int i = 1; i <= loclines.Size(); i++)
		{
			(*testout) << "line ";
			for (int j = 1; j <= 2; j++)
			{
			int hi = 0;
			if (loclines.Get(i).I(j) >= 1 && loclines.Get(i).I(j) <= pindex.Size())
			{
			  hi = adfront.GetGlobalIndex(pindex.Get(loclines.Get(i).I(j)));
			}

			(*testout) << hi << " ";
			}
			(*testout) << " : " << plainpoints.Get(loclines.Get(i).I1()) << " - " << plainpoints.Get(loclines.Get(i).I2()) << " 3d: " << locpoints.Get(loclines.Get(i).I1()) << " - " << locpoints.Get(loclines.Get(i).I2()) << "\n";
		}



		netgen.GlobalMembers.glrender(1);
		}
	}
	else
	{
		adfront.IncrementClass(lindex.Get(1));

		if (debugparam.haltnosuccess || debugflag)
		{
		Console.Write("Problem with seg ");
		Console.Write(gpi1);
		Console.Write(" - ");
		Console.Write(gpi2);
		Console.Write(", class = ");
		Console.Write(qualclass);
		Console.Write("\n");

		(*testout) << "Problem with seg " << gpi1 << " - " << gpi2 << ", class = " << qualclass << "\n";

		multithread.drawing = 1;
		multithread.testmode = 1;
		multithread.pause = 1;


		/*
		  extern STLGeometry * stlgeometry;
		  stlgeometry->ClearMarkedSegs();
		  for (i = 1; i <= loclines.Size(); i++)
		  {
		  stlgeometry->AddMarkedSeg(locpoints.Get(loclines.Get(i).I1()),
		  locpoints.Get(loclines.Get(i).I2()));
		  }
		*/

		for (int i = 1; i <= loclines.Size(); i++)
		{
			(*testout) << "line ";
			for (int j = 1; j <= 2; j++)
			{
			int hi = 0;
			if (loclines.Get(i).I(j) >= 1 && loclines.Get(i).I(j) <= pindex.Size())
			{
			  hi = adfront.GetGlobalIndex(pindex.Get(loclines.Get(i).I(j)));
			}

			(*testout) << hi << " ";
			}
			(*testout) << " : " << plainpoints.Get(loclines.Get(i).I1()) << " - " << plainpoints.Get(loclines.Get(i).I2()) << " 3d: " << locpoints.Get(loclines.Get(i).I1()) << " - " << locpoints.Get(loclines.Get(i).I2()) << "\n";
		}


		/*
		  cout << "p1gi = " << blgeominfo[0].trignum
		  << ", p2gi = " << blgeominfo[1].trignum << endl;
		*/

		netgen.GlobalMembers.glrender(1);
		}


#if MYGRAPH
		if (silentflag < 3)
		{
		if (testmode || trials % 2 == 0)
		{
			plotsurf.DrawPnL(locpoints, loclines);
			plotsurf.Plot(testmode, testmode);
		}
		}
#endif
	}

	}

	PrintMessage(3, "Surface meshing done");


	adfront.PrintOpenSegments(*testout);

	multithread.task = savetask;


	EndMesh();


	if (!adfront.Empty())
	{
	  return MESHING2_GIVEUP;
	}

	return MESHING2_OK;
  }

//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  DLL_HEADER void Delaunay(Mesh mesh, int domainnr, MeshingParameters mp);
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  DLL_HEADER void BlockFillLocalH(Mesh mesh, MeshingParameters mp);


  ///
  public void AddPoint(Point3d p, PointIndex globind, MultiPointGeomInfo mgi = null, bool pointonsurface = true)
  {
	//(*testout) << "add point " << globind << endl;
	adfront.AddPoint(p, globind, mgi, pointonsurface);
  }

  ///
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  DLL_HEADER void AddBoundaryElement(INDEX i1, INDEX i2, PointGeomInfo gi1, PointGeomInfo gi2);

  ///
  public void SetStartTime(double astarttime)
  {
	starttime = astarttime;
  }

  ///
  public void SetMaxArea(double amaxarea)
  {
	maxarea = amaxarea;
  }

  ///
  protected virtual void StartMesh()
  {
	foundmap.SetSize(rules.Size());
	canuse.SetSize(rules.Size());
	ruleused.SetSize(rules.Size());

	foundmap = 0;
	canuse = 0;
	ruleused = 0;

	// cntelem = 0;
	// trials = 0;
  }

  ///
  protected virtual void EndMesh()
  {
	for (int i = 0; i < ruleused.Size(); i++)
	{
	  (*testout) << setw(4) << ruleused[i] << " times used rule " << rules[i].Name() << "\n";
	}
  }

  ///
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual double CalcLocalH(const Point3d & p, double gh) const
  protected virtual double CalcLocalH(Point3d p, double gh)
  {
	return gh;
  }

  ///

  // should be class variables !!(?)
  // static Vec3d ex, ey;
  // static Point3d globp1;

  protected virtual void DefineTransformation(Point3d p1, Point3d p2, PointGeomInfo geominfo1, PointGeomInfo geominfo2)
  {
	globp1 = p1;
	ex = p2 - p1;
	ex /= ex.Length();
	ey.X() = -ex.Y();
	ey.Y() = ex.X();
	ey.Z() = 0;
  }

  ///
  protected virtual void TransformToPlain(Point3d locpoint, MultiPointGeomInfo geominf, Point2d plainpoint, double h, ref int zone)
  {
	Vec3d p1p = new Vec3d(globp1, locpoint);

	//    p1p = locpoint - globp1;
	p1p /= h;
	plainpoint.X() = p1p * ex;
	plainpoint.Y() = p1p * ey;
	zone = 0;
  }

  /// return 0 .. ok
  /// return >0 .. cannot transform point to true surface
  protected virtual int TransformFromPlain(Point2d plainpoint, ref Point3d locpoint, PointGeomInfo gi, double h)
  {
	Vec3d p1p = new Vec3d();
	gi.trignum = 1;

	p1p = plainpoint.X() * ex + plainpoint.Y() * ey;
	p1p *= h;
	locpoint = globp1 + p1p;
	return 0;
  }

  /// projects to surface
  /// return 0 .. ok
  protected virtual int BelongsToActiveChart(Point3d p, PointGeomInfo gi)
  {
	return 1;
  }

  /// computes geoinfo data for line with respect to
  /// selected chart
  protected virtual int ComputePointGeomInfo(Point3d p, PointGeomInfo gi)
  {
	gi.trignum = 1;
	return 0;
  }

  /// Tries to select unique geominfo on active chart
  /// return 0: success
  /// return 1: failed
  protected virtual int ChooseChartPointGeomInfo(MultiPointGeomInfo mpgi, ref PointGeomInfo pgi)
  {
	pgi = mpgi.GetPGI(1);
	return 0;
  }



  /*
    tests, whether endpoint (= 1 or 2) of line segment p1-p2
    is inside of the selected chart. The endpoint must be on the
    chart
   */
  protected virtual int IsLineVertexOnChart(Point3d p1, Point3d p2, int endpoint, PointGeomInfo geominfo)
  {
	return 1;
  }

  /*
    get (projected) boundary of current chart
   */
//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual void GetChartBoundary(Array<Point2d> & points, Array<Point3d> & points3d, Array<INDEX_2> & lines, double h) const
  protected virtual void GetChartBoundary(Array<Point2d> points, Array<Point3d> points3d, Array<INDEX_2> lines, double h)
  {
	points.SetSize(0);
	points3d.SetSize(0);
	lines.SetSize(0);
  }

//C++ TO C# CONVERTER WARNING: 'const' methods are not available in C#:
//ORIGINAL LINE: virtual double Area() const
  protected virtual double Area()
  {
	return -1;
  }


/** Applies 2D rules.
 Tests all 2D rules */
//C++ TO C# CONVERTER TODO TASK: The implementation of the following method could not be found:
//  int ApplyRules(Array<Point2d> lpoints, Array<int> legalpoints, int maxlegalpoint, Array<INDEX_2> llines, int maxlegelline, Array<Element2d> elements, Array<INDEX> dellines, int tolerance, MeshingParameters mp);


}























// #define OPENGL
#if OPENGLxx

/* *********************** Draw Surface Meshing **************** */


#define GL_GLEXT_PROTOTYPES
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
#define GLX_GLXEXT_PROTOTYPES
#define GL_CLAMP_TO_EDGE
#define GL_ARRAY_BUFFER
#define GL_ELEMENT_ARRAY_BUFFER
#define GL_STATIC_DRAW
#define USE_BUFFERS
#define STLBASE



#else
#endif
