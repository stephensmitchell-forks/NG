namespace netgen
{

	public class Refinement
	{
//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void Bisect(Mesh mesh, BisectionOptions opt, Array<double> quality_loss)
		  {
			PrintMessage(1, "Mesh bisection");
			PushStatus("Mesh bisection");
        
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static int timer = NgProfiler::CreateTimer("Bisect");
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static int timer1 = NgProfiler::CreateTimer("Bisect 1");
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static int timer1a = NgProfiler::CreateTimer("Bisect 1a");
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static int timer1b = NgProfiler::CreateTimer("Bisect 1b");
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static int timer2 = NgProfiler::CreateTimer("Bisect 2");
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static int timer2a = NgProfiler::CreateTimer("Bisect 2a");
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static int timer2b = NgProfiler::CreateTimer("Bisect 2b");
			// static int timer2c = NgProfiler::CreateTimer ("Bisect 2c");
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static int timer3 = NgProfiler::CreateTimer("Bisect 3");
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static int timer3a = NgProfiler::CreateTimer("Bisect 3a");
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static int timer3b = NgProfiler::CreateTimer("Bisect 3b");
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static int timer_bisecttet = NgProfiler::CreateTimer("Bisect tets");
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static int timer_bisecttrig = NgProfiler::CreateTimer("Bisect trigs");
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static int timer_bisectsegms = NgProfiler::CreateTimer("Bisect segms");
			NgProfiler.RegionTimer reg1 = new NgProfiler.RegionTimer(Bisect_timer);
        
			opt.tracer("Bisect", false);
        
			NgProfiler.StartTimer(Bisect_timer1);
			NgProfiler.StartTimer(Bisect_timer1a);
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static int localizetimer = NgProfiler::CreateTimer("localize edgepoints");
			NgProfiler.RegionTimer loct = new NgProfiler.RegionTimer(Bisect_localizetimer);
			LocalizeEdgePoints(mesh);
			loct = null;
        
			Array< Array<int,PointIndex.BASE> > idmaps = new Array< Array<int,PointIndex.BASE> >();
			for (int i = 1; i <= mesh.GetIdentifications().GetMaxNr(); i++)
			{
			if (mesh.GetIdentifications().GetType(i) == Identifications.PERIODIC)
			{
				idmaps.Append(new Array<int,PointIndex.BASE>());
				mesh.GetIdentifications().GetMap(i, idmaps.Last(), true);
			}
			}
        
        
			string refelementinfofileread = "";
			string refelementinfofilewrite = "";
        
			if (opt.refinementfilename)
			{
			ifstream inf = new ifstream(opt.refinementfilename);
			string st;
			inf >> st;
			if (st == "refinementinfo")
			{
				while (inf != null)
				{
				while (inf != null && st != "markedelementsfile")
				{
				  inf >> st;
				}
        
				if (inf != null)
				{
				  inf >> st;
				}
        
				if (st == "read" && inf != null)
				{
				  ReadEnclString(inf, refelementinfofileread, '\"');
				}
				else if (st == "write" && inf)
				{
				  ReadEnclString(inf, refelementinfofilewrite, '\"');
				}
				}
			}
			inf.close();
			}
        
        
        
			if (mesh.mglevels == 1 || idmaps.Size() > 0)
			{
			  BisectTetsCopyMesh(mesh, null, opt, idmaps, refelementinfofileread);
			}
        
        
			mesh.ComputeNVertices();
        
			int np = mesh.GetNV();
			mesh.SetNP(np);
        
			// int ne = mesh.GetNE();
			// int nse = mesh.GetNSE();
			// int i, j, l;
        
			// int initnp = np;
			//  int maxsteps = 3;
        
			mesh.mglevels++;
        
			/*
			  if (opt.refinementfilename || opt.usemarkedelements)
			  maxsteps = 3;
			*/
        
        
			if (opt.refine_p)
			{
			int ne = mesh.GetNE();
			int nse = mesh.GetNSE();
			int ox;
			int oy;
			int oz;
			for (ElementIndex ei = 0; ei < ne; ei++)
			{
			  if (mesh[ei].TestRefinementFlag())
			  {
				  mesh[ei].GetOrder(ox,oy,oz);
				  mesh[ei].SetOrder(ox + 1,oy + 1,oz + 1);
				  if (mesh[ei].TestStrongRefinementFlag())
				  {
				mesh[ei].SetOrder(ox + 2,oy + 2,oz + 2);
				  }
			  }
			}
			for (SurfaceElementIndex sei = 0; sei < nse; sei++)
			{
			  if (mesh[sei].TestRefinementFlag())
			  {
				  mesh[sei].GetOrder(ox,oy);
				  mesh[sei].SetOrder(ox + 1,oy + 1);
				  if (mesh[sei].TestStrongRefinementFlag())
				  {
				mesh[sei].SetOrder(ox + 2,oy + 2);
				  }
			  }
			}
        
		#if ! SABINE //Nachbarelemente mit ordx,ordy,ordz
        
			  Array<int,PointIndex.BASE> v_order = new Array<int,PointIndex.BASE>(mesh.GetNP());
			  v_order = 0;
        
			  for (ElementIndex ei = 0; ei < ne; ei++)
			  {
					for (int j = 0; j < mesh[ei].GetNP(); j++)
					{
					  if (mesh[ei].GetOrder() > v_order[mesh[ei][j]])
					  {
						v_order[mesh[ei][j]] = mesh[ei].GetOrder();
					  }
					}
			  }
        
			  for (SurfaceElementIndex sei = 0; sei < nse; sei++)
			  {
					for (int j = 0; j < mesh[sei].GetNP(); j++)
					{
					  if (mesh[sei].GetOrder() > v_order[mesh[sei][j]])
					  {
						v_order[mesh[sei][j]] = mesh[sei].GetOrder();
					  }
					}
			  }
        
			  for (ElementIndex ei = 0; ei < ne; ei++)
			  {
					for (int j = 0; j < mesh[ei].GetNP(); j++)
					{
					  if (mesh[ei].GetOrder() < v_order[mesh[ei][j]] - 1)
					  {
						mesh[ei].SetOrder(v_order[mesh[ei][j]] - 1);
					  }
					}
			  }
        
			  for (SurfaceElementIndex sei = 0; sei < nse; sei++)
			  {
					for (int j = 0; j < mesh[sei].GetNP(); j++)
					{
					  if (mesh[sei].GetOrder() < v_order[mesh[sei][j]] - 1)
					  {
						mesh[sei].SetOrder(v_order[mesh[sei][j]] - 1);
					  }
					}
			  }
        
		#endif
        
				  PopStatus();
				  return;
			}
        
        
        
			// INDEX_2_HASHTABLE<int> cutedges(10 + 5 * (mtets.Size()+mprisms.Size()+mtris.Size()+mquads.Size()));
			INDEX_2_CLOSED_HASHTABLE<PointIndex> cutedges = new INDEX_2_CLOSED_HASHTABLE<PointIndex>(10 + 9 * (mtets.Size() + mprisms.Size() + mtris.Size() + mquads.Size()));
        
			bool noprojection = false;
			NgProfiler.StopTimer(Bisect_timer1a);
        
			for (int l = 1; l <= 1; l++)
			{
			int marked = 0;
			if (opt.refinementfilename)
			{
				ifstream inf = new ifstream(opt.refinementfilename);
				PrintMessage(3, "load refinementinfo from file ", opt.refinementfilename);
        
				string st;
				inf >> st;
				if (st == "refinementinfo")
				{
				  // new version
				for (int i = 1; i <= mtets.Size(); i++)
				{
				  mtets.Elem(i).marked = 0;
				}
				for (int i = 1; i <= mprisms.Size(); i++)
				{
				  mprisms.Elem(i).marked = 0;
				}
				for (int i = 1; i <= mtris.Size(); i++)
				{
				  mtris.Elem(i).marked = 0;
				}
				for (int i = 1; i <= mquads.Size(); i++)
				{
				  mquads.Elem(i).marked = 0;
				}
				for (int i = 1; i <= mprisms.Size(); i++)
				{
				  mids.Elem(i).marked = 0;
				}
        
				inf >> st;
				while (inf != null)
				{
					if (st[0] == '#')
					{
					inf.ignore(10000,'\n');
					inf >> st;
					}
					else if (st == "markedelementsfile")
					{
					inf >> st;
					ReadEnclString(inf, st, '\"');
					inf >> st;
					}
					else if (st == "noprojection")
					{
					noprojection = true;
					inf >> st;
					}
					else if (st == "refine")
					{
					inf >> st;
					if (st == "elements")
					{
						inf >> st;
						bool isint = true;
						for (int ii = 0; isint && ii < st.Length; ii++)
						{
						  isint = (char.IsDigit(st[ii]) != 0);
						}
        
						while (inf != null && isint)
						{
						mtets.Elem(Convert.ToInt32(st)).marked = 3;
						marked = 1;
        
						inf >> st;
						isint = true;
						for (int ii = 0; isint && ii < st.Length; ii++)
						{
						  isint = (char.IsDigit(st[ii]) != 0);
						}
						}
					}
					else if (st == "orthobrick")
					{
						double[] bounds = new double[6];
						for (int i = 0; i < 6; i++)
						{
						  inf >> bounds[i];
						}
        
						int cnt = 0;
        
						for (ElementIndex ei = 0; ei < mesh.GetNE(); ei++)
						{
						Element el = mesh[ei];
        
						//
						Point < 3> center(0,0,0);
						for (int i = 0; i < el.GetNP(); i++)
						{
							MeshPoint point = mesh[el[i]];
							center(0) += point(0);
							center(1) += point(1);
							center(2) += point(2);
						}
						for (int i = 0; i < 3; i++)
						{
						  center(i) *= 1.0 / (double)el.GetNP();
						}
						if (bounds[0] <= center(0) && center(0) <= bounds[3] && bounds[1] <= center(1) && center(1) <= bounds[4] && bounds[2] <= center(2) && center(2) <= bounds[5])
						{
							mtets[ei].marked = 3;
							cnt++;
						}
        
        
		// 				bool contained = false;
		// 				for(int i=0; !contained && i<el.GetNP(); i++)
		// 				  {
		// 				    const MeshPoint & point = mesh[el[i]];
		// 				    contained = (bounds[0] <= point.X() && point.X() <= bounds[3] &&
		// 						 bounds[1] <= point.Y() && point.Y() <= bounds[4] &&
		// 						 bounds[2] <= point.Z() && point.Z() <= bounds[5]);
		// 				  }
		// 				if(contained)
		// 				  {
		// 				    mtets[ei].marked = 3;
		// 				    cnt++;
		// 				  }
						}
        
        
						ostringstream strstr = new ostringstream();
						strstr.precision(2);
		//C++ TO C# CONVERTER TODO TASK: Statements that are interrupted by preprocessor statements are not converted by C++ to C# Converter:
						strstr << "marked " << (float)cnt / (float)mesh.GetNE() * 100.
		#if WIN32
		//C++ TO C# CONVERTER TODO TASK: Statements that are interrupted by preprocessor statements are not converted by C++ to C# Converter:
						   << "%%"
		#else
		//C++ TO C# CONVERTER TODO TASK: Statements that are interrupted by preprocessor statements are not converted by C++ to C# Converter:
						   << "%"
		#endif
						   <<" of the elements";
						PrintMessage(4, strstr.str());
        
						if (cnt > 0)
						{
						  marked = 1;
						}
        
        
						inf >> st;
					}
					else
					{
						throw new Exception("something wrong with refinementinfo file");
					}
					}
				}
				}
				else
				{
				inf.close();
				inf.open(opt.refinementfilename);
        
				char ch;
				for (int i = 1; i <= mtets.Size(); i++)
				{
					inf >> ch;
					if (inf == null)
					{
					  throw new Exception("something wrong with refinementinfo file (old format)");
					}
					mtets.Elem(i).marked = (ch == '1');
				}
				marked = 1;
				}
				inf.close();
			}
        
			else if (opt.usemarkedelements)
			{
				int cntm = 0;
        
				// all in one !
				if (mprisms.Size())
				{
				int cnttet = 0;
				int cntprism = 0;
				for (int i = 1; i <= mesh.GetNE(); i++)
				{
					if (mesh.VolumeElement(i).GetType() == ELEMENT_TYPE.TET || mesh.VolumeElement(i).GetType() == ELEMENT_TYPE.TET10)
					{
					cnttet++;
					mtets.Elem(cnttet).marked = 3 * mesh.VolumeElement(i).TestRefinementFlag();
					if (mtets.Elem(cnttet).marked)
					{
					  cntm++;
					}
					}
					else
					{
					cntprism++;
					mprisms.Elem(cntprism).marked = 2 * mesh.VolumeElement(i).TestRefinementFlag();
					if (mprisms.Elem(cntprism).marked)
					{
					  cntm++;
					}
					}
        
				}
				}
				else
				{
				  for (int i = 1; i <= mtets.Size(); i++)
				  {
				  mtets.Elem(i).marked = 3 * mesh.VolumeElement(i).TestRefinementFlag();
				  if (mtets.Elem(i).marked)
				  {
					cntm++;
				  }
				  }
				}
        
				// (*testout) << "mtets = " << mtets << endl;
        
				/*
				  for (i = 1; i <= mtris.Size(); i++)
				  mtris.Elem(i).marked = 0;
				  for (i = 1; i <= mquads.Size(); i++)
				  mquads.Elem(i).marked = 0;
				*/
        
				if (printmessage_importance > 0)
				{
				ostringstream str = new ostringstream();
				str << "marked elements: " << cntm;
				PrintMessage(4, str.str());
				}
        
				int cnttrig = 0;
				int cntquad = 0;
				for (int i = 1; i <= mesh.GetNSE(); i++)
				{
				if (mesh.SurfaceElement(i).GetType() == ELEMENT_TYPE.TRIG || mesh.SurfaceElement(i).GetType() == ELEMENT_TYPE.TRIG6)
				{
					cnttrig++;
					mtris.Elem(cnttrig).marked = mesh.SurfaceElement(i).TestRefinementFlag() ? 2 : 0;
					// mtris.Elem(cnttrig).marked = 0;
					if (mtris.Elem(cnttrig).marked)
					{
					  cntm++;
					}
				}
				else
				{
					cntquad++;
							// 2d: marked=2, 3d prisms: marked=1
					mquads.Elem(cntquad).marked = mesh.SurfaceElement(i).TestRefinementFlag() ? 4 - mesh.GetDimension() : 0;
					// mquads.Elem(cntquad).marked = 0;
					if (mquads.Elem(cntquad).marked)
					{
					  cntm++;
					}
				}
				}
        
					  if (printmessage_importance > 0)
					  {
				  ostringstream str = new ostringstream();
				  str << "with surface-elements: " << cntm;
				  PrintMessage(4, str.str());
					  }
        
					// he/sz: das wird oben schon richtig gemacht.
					// hier sind die quads vergessen!
					/*
				if (mesh.GetDimension() == 2)
				  {
				cntm = 0;
				for (i = 1; i <= mtris.Size(); i++)
				  {
					mtris.Elem(i).marked =
					  2 * mesh.SurfaceElement(i).TestRefinementFlag();
					//		  mtris.Elem(i).marked = 2;
					if (mtris.Elem(i).marked)
					  cntm++;
				  }
		
				if (!cntm)
				  {
					for (i = 1; i <= mtris.Size(); i++)
					  {
					mtris.Elem(i).marked = 2;
					cntm++;
					  }
				  }
				cout << "trigs: " << mtris.Size() << " ";
				cout << "marked: " << cntm << endl;
				  }
					*/ 
        
				marked = (cntm > 0);
			}
			else
			{
				marked = BTMarkTets(mtets, mprisms, mesh);
			}
        
			if (marked == 0)
			{
				break;
			}
        
        
			//(*testout) << "mtets " << mtets << endl;
        
			if (opt.refine_p)
			{
				PrintMessage(3, "refine p");
        
				for (int i = 1; i <= mtets.Size(); i++)
				{
				  mtets.Elem(i).incorder = mtets.Elem(i).marked ? 1 : 0;
				}
        
				for (int i = 1; i <= mtets.Size(); i++)
				{
				  if (mtets.Elem(i).incorder)
				  {
				mtets.Elem(i).marked = 0;
				  }
				}
        
        
				for (int i = 1; i <= mprisms.Size(); i++)
				{
				  mprisms.Elem(i).incorder = mprisms.Elem(i).marked ? 1 : 0;
				}
        
				for (int i = 1; i <= mprisms.Size(); i++)
				{
				  if (mprisms.Elem(i).incorder)
				  {
				mprisms.Elem(i).marked = 0;
				  }
				}
        
        
				for (int i = 1; i <= mtris.Size(); i++)
				{
				  mtris.Elem(i).incorder = mtris.Elem(i).marked ? 1 : 0;
				}
        
				for (int i = 1; i <= mtris.Size(); i++)
				{
				if (mtris.Elem(i).incorder)
				{
				  mtris.Elem(i).marked = 0;
				}
				}
			}
        
			if (opt.refine_hp)
			{
				PrintMessage(3, "refine hp");
				BitArray singv = new BitArray(np);
				singv.Clear();
        
				if (mesh.GetDimension() == 3)
				{
				for (int i = 1; i <= mesh.GetNSeg(); i++)
				{
					Segment seg = mesh.LineSegment(i);
					singv.Set(new netgen.Segment(seg[0]));
					singv.Set(new netgen.Segment(seg[1]));
				}
				/*
				  for ( i=1; i<= mesh.GetNSE(); i++)
				  {
				  const Element2d & sel = mesh.SurfaceElement(i);
				  for(int j=1; j<=sel.GetNP(); j++)
				  singv.Set(sel.PNum(j));
				  }
				*/
				}
				else
				{
				// vertices with 2 different bnds
				Array<int> bndind = new Array<int>(np);
				bndind = 0;
				for (int i = 1; i <= mesh.GetNSeg(); i++)
				{
					Segment seg = mesh.LineSegment(i);
					for (int j = 0; j < 2; j++)
					{
					int pi = (j == 0) ? seg[0] : seg[1];
					if (bndind.Elem(pi) == 0)
					{
					  bndind.Elem(pi) = seg.edgenr;
					}
					else if (bndind.Elem(pi) != seg.edgenr)
					{
					  singv.Set(pi);
					}
					}
				}
				}
        
        
        
				for (int i = 1; i <= mtets.Size(); i++)
				{
				  mtets.Elem(i).incorder = 1;
				}
				for (int i = 1; i <= mtets.Size(); i++)
				{
				if (!mtets.Elem(i).marked)
				{
				  mtets.Elem(i).incorder = 0;
				}
				for (int j = 0; j < 4; j++)
				{
				  if (singv.Test(mtets.Elem(i).pnums[j]))
				  {
					mtets.Elem(i).incorder = 0;
				  }
				}
				}
				for (int i = 1; i <= mtets.Size(); i++)
				{
				  if (mtets.Elem(i).incorder)
				  {
				mtets.Elem(i).marked = 0;
				  }
				}
        
        
				for (int i = 1; i <= mprisms.Size(); i++)
				{
				  mprisms.Elem(i).incorder = 1;
				}
				for (int i = 1; i <= mprisms.Size(); i++)
				{
				if (!mprisms.Elem(i).marked)
				{
				  mprisms.Elem(i).incorder = 0;
				}
				for (int j = 0; j < 6; j++)
				{
				  if (singv.Test(mprisms.Elem(i).pnums[j]))
				  {
					mprisms.Elem(i).incorder = 0;
				  }
				}
				}
				for (int i = 1; i <= mprisms.Size(); i++)
				{
				  if (mprisms.Elem(i).incorder)
				  {
				mprisms.Elem(i).marked = 0;
				  }
				}
        
        
				for (int i = 1; i <= mtris.Size(); i++)
				{
				  mtris.Elem(i).incorder = 1;
				}
				for (int i = 1; i <= mtris.Size(); i++)
				{
				if (!mtris.Elem(i).marked)
				{
				  mtris.Elem(i).incorder = 0;
				}
				for (int j = 0; j < 3; j++)
				{
				  if (singv.Test(mtris.Elem(i).pnums[j]))
				  {
					mtris.Elem(i).incorder = 0;
				  }
				}
				}
				for (int i = 1; i <= mtris.Size(); i++)
				{
				if (mtris.Elem(i).incorder)
				{
				  mtris.Elem(i).marked = 0;
				}
				}
			}
        
        
        
			int hangingvol;
			int hangingsurf;
			int hangingedge;
        
			//cout << "write?" << endl;
			//string yn;
			//cin >> yn;
        
        
				(*testout) << "refine volume elements" << "\n";
			do
			{
				// refine volume elements
					NgProfiler.StartTimer(Bisect_timer_bisecttet);
					opt.tracer("bisecttet", false);
				int nel = mtets.Size();
				for (int i = 1; i <= nel; i++)
				{
				  if (mtets.Elem(i).marked)
				  {
				  MarkedTet oldtet = mtets.Get(i);
        
				  INDEX_2 edge = new INDEX_2(oldtet.pnums[oldtet.tetedge1], oldtet.pnums[oldtet.tetedge2]);
				  edge.Sort();
        
				  PointIndex newp = new PointIndex();
				  if (cutedges.Used(edge))
				  {
					  newp = cutedges.Get(edge);
				  }
				  else
				  {
					  Point < 3> npt = Center(new mesh.Point(edge.I1()), new mesh.Point(edge.I2()));
							  newp.CopyFrom(mesh.AddPoint(npt));
					  cutedges.Set(edge, newp);
				  }
        
				  MarkedTet newtet1 = new MarkedTet();
				  MarkedTet newtet2 = new MarkedTet();
				  BTBisectTet(oldtet, new netgen.PointIndex(newp), newtet1, newtet2);
        
				  mtets.Elem(i) = newtet1;
				  mtets.Append(newtet2);
        
				  mesh.mlparentelement.Append(i);
				  }
				}
					NgProfiler.StopTimer(Bisect_timer_bisecttet);
					opt.tracer("bisecttet", true);
				int npr = mprisms.Size();
				for (int i = 1; i <= npr; i++)
				{
				  if (mprisms.Elem(i).marked)
				  {
				  MarkedPrism oldprism = new MarkedPrism();
				  MarkedPrism newprism1 = new MarkedPrism();
				  MarkedPrism newprism2 = new MarkedPrism();
				  PointIndex newp1 = new PointIndex();
				  PointIndex newp2 = new PointIndex();
        
				  oldprism = mprisms.Get(i);
				  int pi1 = 0;
				  if (pi1 == oldprism.markededge)
				  {
					pi1++;
				  }
				  int pi2 = 3 - pi1 - oldprism.markededge;
        
				  INDEX_2 edge1 = new INDEX_2(oldprism.pnums[pi1], oldprism.pnums[pi2]);
				  INDEX_2 edge2 = new INDEX_2(oldprism.pnums[pi1 + 3], oldprism.pnums[pi2 + 3]);
				  edge1.Sort();
				  edge2.Sort();
        
				  if (cutedges.Used(edge1))
				  {
					newp1 = cutedges.Get(edge1);
				  }
				  else
				  {
					  Point < 3> npt = Center(new mesh.Point(edge1.I1()), new mesh.Point(edge1.I2()));
							  newp1.CopyFrom(mesh.AddPoint(npt));
					  cutedges.Set(edge1, newp1);
				  }
				  if (cutedges.Used(edge2))
				  {
					newp2 = cutedges.Get(edge2);
				  }
				  else
				  {
					  Point < 3> npt = Center(new mesh.Point(edge2.I1()), new mesh.Point(edge2.I2()));
							  newp2.CopyFrom(mesh.AddPoint(npt));
					  cutedges.Set(edge2, newp2);
				  }
        
        
				  BTBisectPrism(oldprism, new netgen.PointIndex(newp1), new netgen.PointIndex(newp2), newprism1, newprism2);
				  //if(yn == "y")
				  //  (*testout) << "bisected prism " << oldprism << "and got " << newprism1 << "and " << newprism2 << endl;
				  mprisms.Elem(i) = newprism1;
				  mprisms.Append(newprism2);
				  }
				}
        
				int nid = mids.Size();
				for (int i = 1; i <= nid; i++)
				{
				  if (mids.Elem(i).marked)
				  {
				  MarkedIdentification oldid = new MarkedIdentification();
				  MarkedIdentification newid1 = new MarkedIdentification();
				  MarkedIdentification newid2 = new MarkedIdentification();
				  Array<PointIndex> newp = new Array<PointIndex>();
        
				  oldid = mids.Get(i);
        
				  Array<INDEX_2> edges = new Array<INDEX_2>();
				  edges.Append(new INDEX_2(oldid.pnums[oldid.markededge], oldid.pnums[(oldid.markededge+1) % oldid.np]));
				  edges.Append(new INDEX_2(oldid.pnums[oldid.markededge + oldid.np], oldid.pnums[(oldid.markededge+1) % oldid.np + oldid.np]));
        
				  if (oldid.np == 4)
				  {
					  edges.Append(new INDEX_2(oldid.pnums[(oldid.markededge+2) % oldid.np], oldid.pnums[(oldid.markededge+3) % oldid.np]));
					  edges.Append(new INDEX_2(oldid.pnums[(oldid.markededge+2) % oldid.np + oldid.np], oldid.pnums[(oldid.markededge+3) % oldid.np + oldid.np]));
				  }
				  for (int j = 0; j < edges.Size(); j++)
				  {
					  edges[j].Sort();
        
					  if (cutedges.Used(edges[j]))
					  {
					newp.Append(cutedges.Get(edges[j]));
					  }
					  else
					  {
					  Point < 3> npt = Center(new mesh.Point(edges[j].I1()), new mesh.Point(edges[j].I2()));
					  newp.Append(mesh.AddPoint(npt));
					  cutedges.Set(edges[j], newp[j]);
					  }
				  }
        
				  BTBisectIdentification(oldid, newp, newid1, newid2);
				  mids.Elem(i) = newid1;
				  mids.Append(newid2);
				  }
				}
        
        
				//IdentifyCutEdges(mesh, cutedges);
        
					opt.tracer("mark elements", false);
        
				hangingvol = MarkHangingTets(mtets, cutedges, opt.task_manager) + MarkHangingPrisms(mprisms, cutedges) + MarkHangingIdentifications(mids, cutedges);
        
					opt.tracer("mark elements", true);
        
				uint nsel = mtris.Size();
					NgProfiler.StartTimer(Bisect_timer_bisecttrig);
					opt.tracer("Bisect trigs", false);
				for (uint i = 0; i < nsel; i++)
				{
				  if (mtris[i].marked)
				  {
				  MarkedTri newtri1 = new MarkedTri();
				  MarkedTri newtri2 = new MarkedTri();
				  PointIndex newp = new PointIndex();
        
				  MarkedTri oldtri = mtris[i];
				  PointIndex oldpi1 = oldtri.pnums[(oldtri.markededge+1) % 3];
				  PointIndex oldpi2 = oldtri.pnums[(oldtri.markededge+2) % 3];
				  INDEX_2 edge = new INDEX_2(oldpi1, oldpi2);
				  edge.Sort();
        
				  if (cutedges.Used(edge))
				  {
					  newp = cutedges.Get(edge);
				  }
				  else
				  {
					  Point < 3> npt = Center(new mesh.Point(edge.I1()), new mesh.Point(edge.I2()));
					  newp.CopyFrom(mesh.AddPoint(npt));
							  cutedges.Set(edge, newp);
				  }
        
				  int si = mesh.GetFaceDescriptor(oldtri.surfid).SurfNr();
				  PointGeomInfo npgi = new PointGeomInfo();
        
						  if (mesh[newp].Type() != POINTTYPE.EDGEPOINT)
						  {
							PointBetween(new mesh.Point(oldpi1), new mesh.Point(oldpi2), 0.5, si, oldtri.pgeominfo[(oldtri.markededge+1) % 3], oldtri.pgeominfo[(oldtri.markededge+2) % 3], new mesh.Point(newp), npgi);
						  }
        
				  BTBisectTri(oldtri, new netgen.PointIndex(newp), npgi, newtri1, newtri2);
        
				  mtris[i] = newtri1;
				  mtris.Append(newtri2);
				  mesh.mlparentsurfaceelement.Append(i + 1);
				  }
				}
        
					NgProfiler.StopTimer(Bisect_timer_bisecttrig);
					opt.tracer("Bisect trigs", true);
        
				int nquad = mquads.Size();
				for (int i = 1; i <= nquad; i++)
				{
				  if (mquads.Elem(i).marked)
				  {
				  MarkedQuad oldquad = new MarkedQuad();
				  MarkedQuad newquad1 = new MarkedQuad();
				  MarkedQuad newquad2 = new MarkedQuad();
				  PointIndex newp1 = new PointIndex();
				  PointIndex newp2 = new PointIndex();
        
				  oldquad = mquads.Get(i);
						  /*
				  INDEX_2 edge1(oldquad.pnums[0],
						oldquad.pnums[1]);
				  INDEX_2 edge2(oldquad.pnums[2],
						oldquad.pnums[3]);
						  */
						  INDEX_2 edge1 = new INDEX_2();
						  INDEX_2 edge2 = new INDEX_2();
						  PointGeomInfo pgi11 = new PointGeomInfo();
						  PointGeomInfo pgi12 = new PointGeomInfo();
						  PointGeomInfo pgi21 = new PointGeomInfo();
						  PointGeomInfo pgi22 = new PointGeomInfo();
						  if (oldquad.markededge == 0 || oldquad.markededge == 2)
						  {
							edge1.I1() = oldquad.pnums[0];
							pgi11.CopyFrom(oldquad.pgeominfo[0]);
							edge1.I2() = oldquad.pnums[1];
							pgi12.CopyFrom(oldquad.pgeominfo[1]);
							edge2.I1() = oldquad.pnums[2];
							pgi21.CopyFrom(oldquad.pgeominfo[2]);
							edge2.I2() = oldquad.pnums[3];
							pgi22.CopyFrom(oldquad.pgeominfo[3]);
						  }
						  else // 3 || 1
						  {
							edge1.I1() = oldquad.pnums[0];
							pgi11.CopyFrom(oldquad.pgeominfo[0]);
							edge1.I2() = oldquad.pnums[2];
							pgi12.CopyFrom(oldquad.pgeominfo[2]);
							edge2.I1() = oldquad.pnums[1];
							pgi21.CopyFrom(oldquad.pgeominfo[1]);
							edge2.I2() = oldquad.pnums[3];
							pgi22.CopyFrom(oldquad.pgeominfo[3]);
						  }
        
						  edge1.Sort();
				  edge2.Sort();
        
				  if (cutedges.Used(edge1))
				  {
					  newp1 = cutedges.Get(edge1);
				  }
				  else
				  {
					  Point < 3> np1 = Center(new mesh.Point(edge1.I1()), new mesh.Point(edge1.I2()));
							  newp1.CopyFrom(mesh.AddPoint(np1));
					  cutedges.Set(edge1, newp1);
				  }
        
				  if (cutedges.Used(edge2))
				  {
					  newp2 = cutedges.Get(edge2);
				  }
				  else
				  {
					  Point < 3> np2 = Center(new mesh.Point(edge2.I1()), new mesh.Point(edge2.I2()));
					  newp2.CopyFrom(mesh.AddPoint(np2));
					  cutedges.Set(edge2, newp2);
				  }
        
				  PointGeomInfo npgi1 = new PointGeomInfo();
				  PointGeomInfo npgi2 = new PointGeomInfo();
        
				  int si = mesh.GetFaceDescriptor(oldquad.surfid).SurfNr();
				  //		geom->GetSurface(si)->Project (mesh.Point(newp1));
				  //		geom->GetSurface(si)->Project (mesh.Point(newp2));
        
		//                   (*testout)
		//                   cerr << "project point 1 " << newp1 << " old: " << mesh.Point(newp1);
						  PointBetween(new mesh.Point(edge1.I1()), new mesh.Point(edge1.I2()), 0.5, si, pgi11, pgi12, new mesh.Point(newp1), npgi1);
		// 		  (*testout)
		//                   cerr << " new: " << mesh.Point(newp1) << endl;
        
        
		//                   cerr << "project point 2 " << newp2 << " old: " << mesh.Point(newp2);
						  PointBetween(new mesh.Point(edge2.I1()), new mesh.Point(edge2.I2()), 0.5, si, pgi21, pgi22, new mesh.Point(newp2), npgi2);
		//                   cerr << " new: " << mesh.Point(newp2) << endl;
        
        
				  BTBisectQuad(oldquad, new netgen.PointIndex(newp1), npgi1, new netgen.PointIndex(newp2), npgi2, newquad1, newquad2);
        
				  mquads.Elem(i) = newquad1;
				  mquads.Append(newquad2);
				  }
				}
        
					NgProfiler.StartTimer(Bisect_timer1b);
				hangingsurf = MarkHangingTris(mtris, cutedges, opt.task_manager) + MarkHangingQuads(mquads, cutedges);
        
				hangingedge = mesh.GetDimension() == 3 ? 0 : MarkHangingIdentifications(mids, cutedges);
					NgProfiler.StopTimer(Bisect_timer1b);
					NgProfiler.StartTimer(Bisect_timer_bisectsegms);
				int nseg = mesh.GetNSeg();
				for (int i = 1; i <= nseg; i++)
				{
				Segment seg = mesh.LineSegment(i);
				INDEX_2 edge = new INDEX_2(seg[0], seg[1]);
				edge.Sort();
				if (cutedges.Used(edge))
				{
					hangingedge = 1;
					Segment nseg1 = new Segment(seg);
					Segment nseg2 = new Segment(seg);
        
					int newpi = cutedges.Get(edge);
        
					nseg1[1] = newpi;
					nseg2[0] = newpi;
        
					EdgePointGeomInfo newepgi = new EdgePointGeomInfo();
        
        
		//                     
		//                     cerr << "move edgepoint " << newpi << " from " << mesh.Point(newpi);
					PointBetween(new mesh.Point(seg[0]), new mesh.Point(seg[1]), 0.5, seg.surfnr1, seg.surfnr2, seg.epgeominfo[0], seg.epgeominfo[1], new mesh.Point(newpi), newepgi);
		// 		    cerr << " to " << mesh.Point (newpi) << endl;
        
        
					nseg1.epgeominfo[1] = newepgi;
					nseg2.epgeominfo[0] = newepgi;
        
					mesh.LineSegment(i) = nseg1;
					mesh.AddSegment(nseg2);
				}
				}
        
					NgProfiler.StopTimer(Bisect_timer_bisectsegms);
			} while (hangingvol != 0 || hangingsurf != 0 || hangingedge != 0);
        
			/*
				if (printmessage_importance>0)
			  {
			    ostringstream strstr;
			    strstr << mtets.Size() << " tets" << endl
				   << mtris.Size() << " trigs" << endl;
			    if (mprisms.Size())
			      {
				strstr << mprisms.Size() << " prisms" << endl
				       << mquads.Size() << " quads" << endl;
			      }
			    strstr << mesh.GetNP() << " points";
			    PrintMessage(4,strstr.str());
			  }
			*/
			PrintMessage(4, mtets.Size(), " tets");
			PrintMessage(4, mtris.Size(), " trigs");
			if (mprisms.Size())
			{
				PrintMessage(4, mprisms.Size(), " prisms");
				PrintMessage(4, mquads.Size(), " quads");
			}
			PrintMessage(4, mesh.GetNP(), " points");
			}
        
			NgProfiler.StopTimer(Bisect_timer1);
        
        
			// (*testout) << "mtets = " << mtets << endl;
        
			if (opt.refine_hp)
			{
			//
			Array<int> v_order = new Array<int>(mesh.GetNP());
			v_order = 0;
			if (mesh.GetDimension() == 3)
			{
				for (int i = 1; i <= mtets.Size(); i++)
				{
				  if (mtets.Elem(i).incorder)
				  {
				mtets.Elem(i).order++;
				  }
				}
        
				for (int i = 0; i < mtets.Size(); i++)
				{
				  for (int j = 0; j < 4; j++)
				  {
				if ((int)(mtets[i].order) > v_order.Elem(mtets[i].pnums[j]))
				{
				  v_order.Elem(mtets[i].pnums[j]) = mtets[i].order;
				}
				  }
				}
				for (int i = 0; i < mtets.Size(); i++)
				{
				  for (int j = 0; j < 4; j++)
				  {
				if ((int)(mtets[i].order) < v_order.Elem(mtets[i].pnums[j]) - 1)
				{
				  mtets[i].order = v_order.Elem(mtets[i].pnums[j]) - 1;
				}
				  }
				}
			}
			else
			{
				for (int i = 1; i <= mtris.Size(); i++)
				{
				  if (mtris.Elem(i).incorder)
				  {
				  mtris.Elem(i).order++;
				  }
				}
        
				for (int i = 0; i < mtris.Size(); i++)
				{
				  for (int j = 0; j < 3; j++)
				  {
				if ((int)(mtris[i].order) > v_order.Elem(mtris[i].pnums[j]))
				{
				  v_order.Elem(mtris[i].pnums[j]) = mtris[i].order;
				}
				  }
				}
				for (int i = 0; i < mtris.Size(); i++)
				{
				for (int j = 0; j < 3; j++)
				{
				  if ((int)(mtris[i].order) < v_order.Elem(mtris[i].pnums[j]) - 1)
				  {
					mtris[i].order = v_order.Elem(mtris[i].pnums[j]) - 1;
				  }
				}
				}
			}
			}
        
			NgProfiler.StartTimer(Bisect_timer2);
        
			NgProfiler.StartTimer(Bisect_timer2a);
        
			mtets.SetAllocSize(mtets.Size());
			mprisms.SetAllocSize(mprisms.Size());
			mids.SetAllocSize(mids.Size());
			mtris.SetAllocSize(mtris.Size());
			mquads.SetAllocSize(mquads.Size());
        
			opt.tracer("copy tets", false);
			mesh.ClearVolumeElements();
			mesh.VolumeElements().SetAllocSize(mtets.Size() + mprisms.Size());
			mesh.VolumeElements().SetSize(mtets.Size());
			/*
			for (int i = 1; i <= mtets.Size(); i++)
			  {
			Element el(TET);
			el.SetIndex (mtets.Get(i).matindex);
			for (int j = 1; j <= 4; j++)
			  el.PNum(j) = mtets.Get(i).pnums[j-1];
			el.SetOrder (mtets.Get(i).order);
			mesh.AddVolumeElement (el);
			  }
			*/
			ParallelForRange(opt.task_manager, mtets.Size(), (uint begin, uint end) =>
			{
				 for (uint i = begin; i < end; i++)
				 {
					Element el = new Element(ELEMENT_TYPE.TET);
					var tet = mtets[i];
					el.SetIndex(tet.matindex);
					el.SetOrder(tet.order);
					for (int j = 0; j < 4; j++)
					{
					  el[j] = tet.pnums[j];
					}
					mesh.SetVolumeElement(new ElementIndex((int)i), el);
				 }
			});
        
			opt.tracer("copy tets", true);
        
			for (int i = 1; i <= mprisms.Size(); i++)
			{
			Element el = new Element(ELEMENT_TYPE.PRISM);
			el.SetIndex(mprisms.Get(i).matindex);
			for (int j = 1; j <= 6; j++)
			{
			  el.PNum(j) = mprisms.Get(i).pnums[j - 1];
			}
			el.SetOrder(mprisms.Get(i).order);
        
			// degenerated prism ?
			int[] map1 = {3, 2, 5, 6, 1};
			int[] map2 = {1, 3, 6, 4, 2};
			int[] map3 = {2, 1, 4, 5, 3};
        
        
		//C++ TO C# CONVERTER TODO TASK: C# does not have an equivalent to pointers to value types:
		//ORIGINAL LINE: const int * map = null;
			int map = null;
			int deg1 = 0;
			int deg2 = 0;
			int deg3 = 0;
			// int deg = 0;
			if (el.PNum(1) == el.PNum(4))
			{
				map = map1;
				deg1 = 1;
			}
			if (el.PNum(2) == el.PNum(5))
			{
				map = map2;
				deg2 = 1;
			}
			if (el.PNum(3) == el.PNum(6))
			{
				map = map3;
				deg3 = 1;
			}
        
			switch (deg1 + deg2 + deg3)
			{
			  case 1:
			  {
				  for (int j = 1; j <= 5; j++)
				  {
				el.PNum(j) = mprisms.Get(i).pnums[map[j - 1] - 1];
				  }
        
				  el.SetType(ELEMENT_TYPE.PYRAMID);
				  break;
			  }
			  case 2:
			  {
				  int[] tetmap1 = {1, 2, 3, 4};
				  int[] tetmap2 = {2, 3, 1, 5};
				  int[] tetmap3 = {3, 1, 2, 6};
				  if (deg1 == 0)
				  {
					  map = tetmap1;
				  }
				  if (deg2 == 0)
				  {
					  map = tetmap2;
				  }
				  if (deg3 == 0)
				  {
					  map = tetmap3;
				  }
				  for (int j = 1; j <= 4; j++)
				  {
				el.PNum(j) = mprisms.Get(i).pnums[map[j - 1] - 1];
				  }
				  /*
				if (!deg1) el.PNum(4) = el.PNum(4);
				if (!deg2) el.PNum(4) = el.PNum(5);
				if (!deg3) el.PNum(4) = el.PNum(6);
				  */
				  el.SetType(ELEMENT_TYPE.TET);
				  break;
			  }
			  default:
				;
				break;
			}
			mesh.AddVolumeElement(el);
			}
        
			mesh.ClearSurfaceElements();
			mesh.SurfaceElements().SetAllocSize(mtris.Size() + mquads.Size());
			NgProfiler.StopTimer(Bisect_timer2a);
			NgProfiler.StartTimer(Bisect_timer2b);
			/*
			for (auto & trig : mtris)
			  {
			Element2d el(TRIG);
			el.SetIndex (trig.surfid);
			el.SetOrder (trig.order);
			for (int j = 0; j < 3; j++)
			  {
				el[j] = trig.pnums[j];
				el.GeomInfoPi(j+1) = trig.pgeominfo[j];
			  }
			mesh.AddSurfaceElement (el);
			  }
			*/
        
			mesh.SurfaceElements().SetSize(mtris.Size());
			// for (size_t i = 0; i < mtris.Size(); i++)
			ParallelForRange(opt.task_manager, mtris.Size(), (uint begin, uint end) =>
			{
				 for (uint i = begin; i < end; i++)
				 {
					Element2d el = new Element2d(ELEMENT_TYPE.TRIG);
					var trig = mtris[i];
					el.SetIndex(trig.surfid);
					el.SetOrder(trig.order);
					for (int j = 0; j < 3; j++)
					{
						el[j] = trig.pnums[j];
						el.GeomInfoPi(j + 1) = trig.pgeominfo[j];
					}
					mesh.SetSurfaceElement(new SurfaceElementIndex((int)i), el);
				 }
			});
			mesh.RebuildSurfaceElementLists();
        
			for (int i = 1; i <= mquads.Size(); i++)
			{
			Element2d el = new Element2d(ELEMENT_TYPE.QUAD);
			el.SetIndex(mquads.Get(i).surfid);
			for (int j = 1; j <= 4; j++)
			{
			  el.PNum(j) = mquads.Get(i).pnums[j - 1];
			}
			Swap(ref el.PNum(3), ref el.PNum(4));
			mesh.AddSurfaceElement(el);
			}
			NgProfiler.StopTimer(Bisect_timer2b);
        
        
			// write multilevel hierarchy to mesh:
			np = mesh.GetNP();
			mesh.mlbetweennodes.SetSize(np);
			if (mesh.mglevels <= 2)
			{
			PrintMessage(4, "RESETTING mlbetweennodes");
			for (int i = 1; i <= np; i++)
			{
				mesh.mlbetweennodes.Elem(i).I1() = 0;
				mesh.mlbetweennodes.Elem(i).I2() = 0;
			}
			}
        
			/*
			  for (i = 1; i <= cutedges.GetNBags(); i++)
			  for (j = 1; j <= cutedges.GetBagSize(i); j++)
			  {
			  INDEX_2 edge;
			  int newpi;
			  cutedges.GetData (i, j, edge, newpi);
			  mesh.mlbetweennodes.Elem(newpi) = edge;
			  }
			*/
        
			BitArray isnewpoint = new BitArray(np);
			isnewpoint.Clear();
        
			for (int i = 0; i < cutedges.Size(); i++)
			{
			  if (cutedges.UsedPos0(i))
			  {
			  INDEX_2 edge = new INDEX_2();
			  PointIndex newpi = new PointIndex();
			  cutedges.GetData0(i, ref edge, ref newpi);
			  isnewpoint.Set(new netgen.PointIndex(newpi));
			  mesh.mlbetweennodes.Elem(newpi) = edge;
			  }
			}
        
        
			/*
			  mesh.PrintMemInfo (cout);
			  cout << "tets ";
			  mtets.PrintMemInfo (cout);
			  cout << "prims ";
			  mprisms.PrintMemInfo (cout);
			  cout << "tris ";
			  mtris.PrintMemInfo (cout);
			  cout << "quads ";
			  mquads.PrintMemInfo (cout);
			  cout << "cutedges ";
			  cutedges.PrintMemInfo (cout);
			*/
        
        
			/*
		
			// find connected nodes (close nodes)
			TABLE<int> conto(np);
			for (i = 1; i <= mprisms.Size(); i++)
			for (j = 1; j <= 6; j++)
			{
			int n1 = mprisms.Get(i).pnums[j-1];
			int n2 = mprisms.Get(i).pnums[(j+2)%6];
			//	    if (n1 != n2)
			{
			int found = 0;
			for (k = 1; k <= conto.EntrySize(n1); k++)
			if (conto.Get(n1, k) == n2)
			{
			found = 1;
			break;
			}
			if (!found)
			conto.Add (n1, n2);
			}
			}
			mesh.connectedtonode.SetSize(np);
			for (i = 1; i <= np; i++)
			mesh.connectedtonode.Elem(i) = 0;
		
		
			//       (*testout) << "connection table: " << endl;
			//       for (i = 1; i <= np; i++)
			//       {
			//       (*testout) << "node " << i << ": ";
			// 	  for (j = 1; j <= conto.EntrySize(i); j++)
			// 	  (*testout) << conto.Get(i, j) << " ";
			// 	  (*testout) << endl;
			// 	}
		
		
			for (i = 1; i <= np; i++)
			if (mesh.connectedtonode.Elem(i) == 0)
			{
			mesh.connectedtonode.Elem(i) = i;
			ConnectToNodeRec (i, i, conto, mesh.connectedtonode);
			}
			*/  
        
			//  mesh.BuildConnectedNodes();
        
        
        
        
			mesh.ComputeNVertices();
        
			opt.tracer("call RebuildSurfElList", false);
			mesh.RebuildSurfaceElementLists();
			opt.tracer("call RebuildSurfElList", true);
        
        
			// update identification tables
			for (int i = 1; i <= mesh.GetIdentifications().GetMaxNr(); i++)
			{
			Array<int,PointIndex.BASE> identmap = new Array<int,PointIndex.BASE>();
        
			mesh.GetIdentifications().GetMap(i, identmap);
        
        
			/*
			  for (j = 1; j <= cutedges.GetNBags(); j++)
			  for (k = 1; k <= cutedges.GetBagSize(j); k++)
			  {
			  INDEX_2 i2;
			  int newpi;
			  cutedges.GetData (j, k, i2, newpi);
			  INDEX_2 oi2(identmap.Get(i2.I1()),
			  identmap.Get(i2.I2()));
			  oi2.Sort();
			  if (cutedges.Used (oi2))
			  {
			  int onewpi = cutedges.Get(oi2);
			  mesh.GetIdentifications().Add (newpi, onewpi, i);
			  }
			  }
			*/
        
			for (int j = 0; j < cutedges.Size(); j++)
			{
			  if (cutedges.UsedPos0(j))
			  {
				  INDEX_2 i2 = new INDEX_2();
				  PointIndex newpi = new PointIndex();
				  cutedges.GetData0(j, ref i2, ref newpi);
				  INDEX_2 oi2 = new INDEX_2(identmap.Get(i2.I1()), identmap.Get(i2.I2()));
				  oi2.Sort();
				  if (cutedges.Used(oi2))
				  {
				  PointIndex onewpi = cutedges.Get(oi2);
				  mesh.GetIdentifications().Add(new netgen.PointIndex(newpi), new netgen.PointIndex(onewpi), i);
				  }
			  }
			}
			}
        
			opt.tracer("Bisect", true);
        
			// Repair works only for tets!
			bool do_repair = mesh.PureTetMesh();
        
			do_repair = false; // JS, March 2009: multigrid crashes
        
			//if(mesh.mglevels == 3)
			//  noprojection = true;
        
			//noprojection = true;
        
			if (noprojection)
			{
			do_repair = false;
			for (int ii = 1; ii <= mesh.GetNP(); ii++)
			{
				if (isnewpoint.Test(ii) && mesh.mlbetweennodes[ii][0] > 0)
				{
				mesh.Point(ii) = Center(new mesh.Point(mesh.mlbetweennodes[ii][0]), new mesh.Point(mesh.mlbetweennodes[ii][1]));
				}
			}
			}
        
        
			// Check/Repair
        
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static bool repaired_once;
			if (mesh.mglevels == 1)
			{
			  Bisect_repaired_once = false;
			}
        
			//mesh.Save("before.vol");
        
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static int reptimer = NgProfiler::CreateTimer("check/repair");
			NgProfiler.RegionTimer regt = new NgProfiler.RegionTimer(null);
			regt = new NgProfiler.RegionTimer(Bisect_reptimer);
        
			Array<ElementIndex> bad_elts = new Array<ElementIndex>();
			Array<double> pure_badness = new Array<double>();
        
			if (do_repair || quality_loss != null)
			{
			pure_badness.SetSize(mesh.GetNP() + 2);
			GetPureBadness(mesh, pure_badness, isnewpoint);
			}
        
        
			if (do_repair) // by Markus W
			{
			const double max_worsening = 1;
        
			const bool uselocalworsening = false;
        
			// bool repaired = false;
        
			Validate(mesh, bad_elts, pure_badness, max_worsening, uselocalworsening);
        
				if (printmessage_importance > 0)
				{
				ostringstream strstr = new ostringstream();
				for (int ii = 0; ii < bad_elts.Size(); ii++)
				{
				  strstr << "bad element " << bad_elts[ii] << "\n";
				}
				PrintMessage(1, strstr.str());
				}
			if (Bisect_repaired_once || bad_elts.Size() > 0)
			{
				clock_t t1 = new clock_t(clock());
        
        
				// update id-maps
				int j = 0;
				for (int i = 1; i <= mesh.GetIdentifications().GetMaxNr(); i++)
				{
				if (mesh.GetIdentifications().GetType(i) == Identifications.PERIODIC)
				{
					mesh.GetIdentifications().GetMap(i, idmaps[j], true);
					j++;
				}
				}
        
        
				// do the repair
				try
				{
				RepairBisection(mesh, bad_elts, isnewpoint, this, pure_badness, max_worsening, uselocalworsening, idmaps);
				// repaired = true;
				Bisect_repaired_once = true;
				}
				catch (Exception ex)
				{
				PrintMessage(1, "Problem: " + ex.What());
				}
        
        
					if (printmessage_importance > 0)
					{
				  ostringstream strstr = new ostringstream();
					  strstr << "Time for Repair: " << (double)(clock() - t1) / (double)DefineConstants.CLOCKS_PER_SEC << "\n" << "bad elements after repair: " << bad_elts << "\n";
				  PrintMessage(1, strstr.str());
					}
        
				if (quality_loss != null)
				{
				  Validate(mesh, bad_elts, pure_badness, 1e100, uselocalworsening, quality_loss);
				}
        
				if (idmaps.Size() == 0)
				{
				  UpdateEdgeMarks(mesh, idmaps);
				}
        
				/*
				if(1==1)
				  UpdateEdgeMarks(mesh,idmaps);
				else
				  mesh.mglevels = 1;
				*/
        
				//mesh.ImproveMesh();
        
			}
			}
			regt = null;
        
        
        
			for (int i = 0; i < idmaps.Size(); i++)
			{
			  idmaps[i] = null;
			}
			idmaps.DeleteAll();
        
			NgProfiler.StopTimer(Bisect_timer2);
			NgProfiler.StartTimer(Bisect_timer3);
        
			NgProfiler.StartTimer(Bisect_timer3a);
			opt.tracer("topology from bisect", false);
			mesh.UpdateTopology(opt.task_manager, opt.tracer);
			opt.tracer("topology from bisect", true);
			NgProfiler.StopTimer(Bisect_timer3a);
        
        
        
			if (refelementinfofilewrite != "")
			{
			PrintMessage(3, "writing marked-elements information to \"", refelementinfofilewrite, "\"");
			ofstream ofst = new ofstream(refelementinfofilewrite);
        
			WriteMarkedElements(ofst);
        
			ofst.close();
			}
        
			NgProfiler.StartTimer(Bisect_timer3b);
			mesh.CalcSurfacesOfNode();
			NgProfiler.StopTimer(Bisect_timer3b);
        
			PrintMessage(1, "Bisection done");
        
			PopStatus();
			NgProfiler.StopTimer(Bisect_timer3);
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void Refine(Mesh mesh)
		  {
		//C++ TO C# CONVERTER TODO TASK: There is no equivalent to 'const_cast' in C#:
			const_cast<Refinement&> (this).Refine(mesh);
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void Refine(Mesh mesh)
		  {
			PrintMessage(3, "Refine mesh");
        
			mesh.SetNextMajorTimeStamp();
        
			if (ntasks > 1 && id == 0)
			{
			  return;
			}
        
        
			// reduce 2nd order
			mesh.ComputeNVertices();
			mesh.SetNP(mesh.GetNV());
        
			if (mesh.mlbetweennodes.Size() < mesh.GetNV())
			{
				mesh.mlbetweennodes.SetSize(mesh.GetNV());
				mesh.mlbetweennodes = new INDEX_2(PointIndex.BASE-1, PointIndex.BASE-1);
			}
        
        
			INDEX_2_HASHTABLE<PointIndex> between = new INDEX_2_HASHTABLE<PointIndex>(mesh.GetNP() + 5);
        
        
			// new version with consistent ordering across sub-domains
        
			Array<INDEX_2> parents = new Array<INDEX_2>();
			for (SegmentIndex si = 0; si < mesh.GetNSeg(); si++)
			{
			Segment el = mesh[si];
			INDEX_2 i2 = INDEX_2.Sort(new netgen.Segment(el[0]), new netgen.Segment(el[1]));
				if (!between.Used(i2))
				{
					between.Set(i2, 0);
					parents.Append(i2);
				}
			}
			for (SurfaceElementIndex sei = 0; sei < mesh.GetNSE(); sei++)
			{
			Element2d el = mesh[sei];
			switch (el.GetType())
			{
			  case ELEMENT_TYPE.TRIG:
			  case ELEMENT_TYPE.TRIG6:
			  {
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//			  static int betw[3][3] = { { 1, 2, 3 }, { 0, 2, 4 }, { 0, 1, 5 } };
					  for (int j = 0; j < 3; j++)
					  {
						  var i2 = PointIndices < 2>.Sort(el[Refine_betw[j][0]],el[Refine_betw[j][1]]);
						  if (!between.Used(i2))
						  {
							  between.Set(i2, 0);
							  parents.Append(i2);
						  }
					  }
					  break;
			  }
				  default:
					throw new Exception("currently refinement for quad-elements is not supported");
			}
			}
			for (ElementIndex ei = 0; ei < mesh.GetNE(); ei++)
			{
			Element el = mesh[ei];
			switch (el.GetType())
			{
			  case ELEMENT_TYPE.TET:
			  case ELEMENT_TYPE.TET10:
			  {
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//		  static int betw[6][3] = { { 1, 2, 5 }, { 1, 3, 6 }, { 1, 4, 7 }, { 2, 3, 8 }, { 2, 4, 9 }, { 3, 4, 10 } };
        
					  for (int j = 0; j < 6; j++)
					  {
						  INDEX_2 i2 = INDEX_2.Sort(el.PNum(Refine_betw[j][0]), el.PNum(Refine_betw[j][1]));
						  if (!between.Used(i2))
						  {
							  between.Set(i2, 0);
							  parents.Append(i2);
						  }
					  }
					  break;
			  }
				  default:
					throw new Exception("currently refinement for non-tet elements is not supported");
			}
			}
        
			PrintMessage(5, "have points");
        
			Array<int> par_nr = new Array<int>(parents.Size());
			for (int i = 0; i < par_nr.Size(); i++)
			{
			  par_nr[i] = i;
			}
			QuickSort(parents, par_nr);
			mesh.mlbetweennodes.SetSize(mesh.GetNV() + parents.Size());
			for (int i = 0; i < parents.Size(); i++)
			{
				between.Set(parents[i], mesh.GetNV() + i + PointIndex.BASE);
				mesh.mlbetweennodes[mesh.GetNV() + i + PointIndex.BASE] = parents[i];
			}
        
			mesh.SetNP(mesh.GetNV() + parents.Size());
			Array<bool, PointIndex.BASE> pointset = new Array<bool, PointIndex.BASE>(mesh.GetNP());
			pointset = false;
        
			PrintMessage(5, "sorting complete");
        
			// refine edges
			Array<EdgePointGeomInfo,PointIndex.BASE> epgi = new Array<EdgePointGeomInfo,PointIndex.BASE>();
        
			int oldns = mesh.GetNSeg();
			for (SegmentIndex si = 0; si < oldns; si++)
			{
		//C++ TO C# CONVERTER TODO TASK: The following line was determined to contain a copy constructor call - this should be verified and a copy constructor should be created:
		//ORIGINAL LINE: const Segment & el = mesh.LineSegment(si);
			Segment el = mesh.LineSegment(new netgen.SegmentIndex(si));
        
			INDEX_2 i2 = INDEX_2.Sort(new netgen.Segment(el[0]), new netgen.Segment(el[1]));
			PointIndex pinew = between.Get(i2);
			EdgePointGeomInfo ngi = new EdgePointGeomInfo();
        
			if (pointset[pinew])
			{
				// pinew = between.Get(i2);
				ngi = epgi[pinew];
			}
			else
			{
					pointset[pinew] = true;
				Point < 3> pnew;
				PointBetween(new mesh.Point(el[0]), new mesh.Point(el[1]), 0.5, el.surfnr1, el.surfnr2, el.epgeominfo[0], el.epgeominfo[1], pnew, ngi);
        
				// pinew = mesh.AddPoint (pnew);
					mesh.Point(new netgen.PointIndex(pinew)) = pnew;
				// between.Set (i2, pinew);
        
				if (pinew >= epgi.Size() + PointIndex.BASE)
				{
				  epgi.SetSize(pinew + 1 - PointIndex.BASE);
				}
				epgi[pinew] = ngi;
			}
        
			Segment ns1 = new Segment(el);
			Segment ns2 = new Segment(el);
			ns1[1] = pinew;
			ns1.epgeominfo[1] = ngi;
			ns2[0] = pinew;
			ns2.epgeominfo[0] = ngi;
        
		//C++ TO C# CONVERTER TODO TASK: The following line was determined to contain a copy constructor call - this should be verified and a copy constructor should be created:
		//ORIGINAL LINE: mesh.LineSegment(si) = ns1;
			mesh.LineSegment(new netgen.SegmentIndex(si)) = ns1;
			mesh.AddSegment(ns2);
			}
        
			PrintMessage(5, "have 1d elements");
        
			// refine surface elements
			Array<PointGeomInfo,PointIndex.BASE> surfgi = new Array<PointGeomInfo,PointIndex.BASE>(8 * mesh.GetNP());
			for (int i = PointIndex.BASE; i < surfgi.Size() + PointIndex.BASE; i++)
			{
			  surfgi[i].trignum = -1;
			}
        
        
			int oldnf = mesh.GetNSE();
			for (SurfaceElementIndex sei = 0; sei < oldnf; sei++)
			{
		//C++ TO C# CONVERTER TODO TASK: The following line was determined to contain a copy constructor call - this should be verified and a copy constructor should be created:
		//ORIGINAL LINE: const Element2d & el = mesh.SurfaceElement(sei);
			Element2d el = mesh.SurfaceElement(new netgen.SurfaceElementIndex(sei));
        
			switch (el.GetType())
			{
			  case ELEMENT_TYPE.TRIG:
			  case ELEMENT_TYPE.TRIG6:
			  {
				  ArrayMem<PointIndex,6> pnums = new ArrayMem<PointIndex,6>(6);
				  ArrayMem<PointGeomInfo,6> pgis = new ArrayMem<PointGeomInfo,6>(6);
        
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//		  static int betw[3][3] = { { 2, 3, 4 }, { 1, 3, 5 }, { 1, 2, 6 } };
        
				  for (int j = 1; j <= 3; j++)
				  {
				  pnums.Elem(j) = el.PNum(j);
				  pgis.Elem(j) = el.GeomInfoPi(j);
				  }
        
				  for (int j = 0; j < 3; j++)
				  {
				  PointIndex pi1 = pnums.Elem(Refine_betw[j][0]);
				  PointIndex pi2 = pnums.Elem(Refine_betw[j][1]);
        
				  INDEX_2 i2 = new INDEX_2(pi1, pi2);
				  i2.Sort();
        
				  Point < 3> pb;
				  PointGeomInfo pgi = new PointGeomInfo();
				  PointBetween(new mesh.Point(pi1), new mesh.Point(pi2), 0.5, mesh.GetFaceDescriptor(el.GetIndex()).SurfNr(), el.GeomInfoPi(Refine_betw[j][0]), el.GeomInfoPi(Refine_betw[j][1]), pb, pgi);
        
        
				  pgis.Elem(4 + j) = pgi;
						  PointIndex pinew = between.Get(i2);
						  pnums.Elem(4 + j) = pinew;
						  if (!pointset[pinew])
						  {
							  pointset[pinew] = true;
							  mesh.Point(new netgen.PointIndex(pinew)) = pb;
						  }
						  /*
				  if (between.Used(i2))
					pnums.Elem(4+j) = between.Get(i2);
				  else
					{
					  pnums.Elem(4+j) = mesh.AddPoint (pb);
					  between.Set (i2, pnums.Get(4+j));
					}
						  */
				  if (surfgi.Size() < pnums.Elem(4 + j))
				  {
					surfgi.SetSize(pnums.Elem(4 + j));
				  }
				  surfgi.Elem(pnums.Elem(4 + j)) = pgis.Elem(4 + j);
				  }
        
        
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//		  static int reftab[4][3] = { { 1, 6, 5 }, { 2, 4, 6 }, { 3, 5, 4 }, { 6, 4, 5 } };
        
				  int ind = el.GetIndex();
				  for (int j = 0; j < 4; j++)
				  {
				  Element2d nel = new Element2d(ELEMENT_TYPE.TRIG);
				  for (int k = 1; k <= 3; k++)
				  {
					  nel.PNum(k) = pnums.Get(Refine_reftab[j][k - 1]);
					  nel.GeomInfoPi(k) = pgis.Get(Refine_reftab[j][k - 1]);
				  }
				  nel.SetIndex(ind);
        
				  if (j == 0)
				  {
		//C++ TO C# CONVERTER TODO TASK: The following line was determined to contain a copy constructor call - this should be verified and a copy constructor should be created:
		//ORIGINAL LINE: mesh.SurfaceElement(sei) = nel;
					mesh.SurfaceElement(new netgen.SurfaceElementIndex(sei)) = nel;
				  }
				  else
				  {
					mesh.AddSurfaceElement(nel);
				  }
				  }
				  break;
			  }
			  case ELEMENT_TYPE.QUAD:
			  case ELEMENT_TYPE.QUAD6:
			  case ELEMENT_TYPE.QUAD8:
			  {
				  ArrayMem<PointIndex,9> pnums = new ArrayMem<PointIndex,9>(9);
				  ArrayMem<PointGeomInfo,9> pgis = new ArrayMem<PointGeomInfo,9>(9);
        
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//		  static int betw[5][3] = { { 1, 2, 5 }, { 2, 3, 6 }, { 3, 4, 7 }, { 1, 4, 8 }, { 5, 7, 9 } };
        
				  for (int j = 1; j <= 4; j++)
				  {
				  pnums.Elem(j) = el.PNum(j);
				  pgis.Elem(j) = el.GeomInfoPi(j);
				  }
        
				  for (int j = 0; j < 5; j++)
				  {
				  int pi1 = pnums.Elem(Refine_betw[j][0]);
				  int pi2 = pnums.Elem(Refine_betw[j][1]);
        
				  INDEX_2 i2 = new INDEX_2(pi1, pi2);
				  i2.Sort();
        
				  if (between.Used(i2))
				  {
					  pnums.Elem(5 + j) = between.Get(i2);
					  pgis.Elem(5 + j) = surfgi.Get(pnums.Elem(4 + j));
				  }
				  else
				  {
					  Point < 3> pb;
					  PointBetween(new mesh.Point(pi1), new mesh.Point(pi2), 0.5, mesh.GetFaceDescriptor(el.GetIndex()).SurfNr(), el.GeomInfoPi(Refine_betw[j][0]), el.GeomInfoPi(Refine_betw[j][1]), pb, pgis.Elem(5 + j));
        
					  pnums.Elem(5 + j) = mesh.AddPoint(pb);
        
					  between.Set(i2, pnums.Get(5 + j));
        
					  if (surfgi.Size() < pnums.Elem(5 + j))
					  {
					surfgi.SetSize(pnums.Elem(5 + j));
					  }
					  surfgi.Elem(pnums.Elem(5 + j)) = pgis.Elem(5 + j);
				  }
				  }
        
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//		  static int reftab[4][4] = { { 1, 5, 9, 8 }, { 5, 2, 6, 9 }, { 8, 9, 7, 4 }, { 9, 6, 3, 7 } };
        
				  int ind = el.GetIndex();
				  for (int j = 0; j < 4; j++)
				  {
				  Element2d nel = new Element2d(ELEMENT_TYPE.QUAD);
				  for (int k = 1; k <= 4; k++)
				  {
					  nel.PNum(k) = pnums.Get(Refine_reftab[j][k - 1]);
					  nel.GeomInfoPi(k) = pgis.Get(Refine_reftab[j][k - 1]);
				  }
				  nel.SetIndex(ind);
        
				  if (j == 0)
				  {
		//C++ TO C# CONVERTER TODO TASK: The following line was determined to contain a copy constructor call - this should be verified and a copy constructor should be created:
		//ORIGINAL LINE: mesh.SurfaceElement(sei) = nel;
					mesh.SurfaceElement(new netgen.SurfaceElementIndex(sei)) = nel;
				  }
				  else
				  {
					mesh.AddSurfaceElement(nel);
				  }
				  }
				  break;
			  }
			  default:
				PrintSysError("Refine: undefined surface element type ", (int)el.GetType());
				break;
			}
			}
        
			PrintMessage(5, "have 2d elements");
			// cout << "id = " << id << ", ne = " << mesh.GetNE() << endl;
			// refine volume elements
			int oldne = mesh.GetNE();
			mesh.VolumeElements().SetAllocSize(8 * oldne);
			for (ElementIndex ei = 0; ei < oldne; ei++)
			{
		//C++ TO C# CONVERTER TODO TASK: The following line was determined to contain a copy constructor call - this should be verified and a copy constructor should be created:
		//ORIGINAL LINE: const Element & el = mesh.VolumeElement(ei);
			Element el = mesh.VolumeElement(new netgen.ElementIndex(ei));
			switch (el.GetType())
			{
			  case ELEMENT_TYPE.TET:
			  case ELEMENT_TYPE.TET10:
			  {
				 ArrayMem<PointIndex,10> pnums = new ArrayMem<PointIndex,10>(10);
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//		 static int betw[6][3] = { { 1, 2, 5 }, { 1, 3, 6 }, { 1, 4, 7 }, { 2, 3, 8 }, { 2, 4, 9 }, { 3, 4, 10 } };
        
				 int elrev = el.flags.reverse;
        
				 for (int j = 1; j <= 4; j++)
				 {
					   pnums.Elem(j) = el.PNum(j);
				 }
				 if (elrev != 0)
				 {
				 swap(pnums.Elem(3), pnums.Elem(4));
				 }
        
				 for (int j = 0; j < 6; j++)
				 {
						 PointIndex pi1 = pnums.Get(Refine_betw[j][0]);
						 PointIndex pi2 = pnums.Get(Refine_betw[j][1]);
						 INDEX_2 i2 = new INDEX_2(pi1, pi2);
						 i2.Sort();
        
				   /*
				   if (between.Used(i2))
				      pnums.Elem(5+j) = between.Get(i2);
				   else
				   {
				  pnums.Elem(5+j) = mesh.AddPoint
				  (Center (mesh.Point(i2.I1()),
					   mesh.Point(i2.I2())));
				  between.Set (i2, pnums.Elem(5+j));
				   }
				   */
				   PointIndex pinew = between.Get(i2);
				   pnums.Elem(j + 5) = pinew;
				   if (!pointset[pinew])
				   {
				   pointset[pinew] = true;
				   mesh.Point(new netgen.PointIndex(pinew)) = Center(new mesh.Point(pi1), new mesh.Point(pi2));
				   }
				 }
        
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//		static int reftab[8][4] = { { 1, 5, 6, 7 }, { 5, 2, 8, 9 }, { 6, 8, 3, 10 }, { 7, 9, 10, 4 }, { 5, 6, 7, 9 }, { 5, 6, 9, 8 }, { 6, 7, 9, 10 }, { 6, 8, 10, 9 } };
			/*
			  { { 1, 5, 6, 7 },
			  { 5, 2, 8, 9 },
			  { 6, 8, 3, 10 },
			  { 7, 9, 10, 4 },
			  { 5, 6, 7, 9 },
			  { 5, 6, 8, 9 },
			  { 6, 7, 9, 10 },
			  { 6, 8, 9, 10 } };
			*/
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	   static bool reverse[8] = { false, false, false, false, false, true, false, true };
        
			   int ind = el.GetIndex();
			   for (int j = 0; j < 8; j++)
			   {
					 Element nel = new Element(ELEMENT_TYPE.TET);
				  for (int k = 1; k <= 4; k++)
				  {
					nel.PNum(k) = pnums.Get(Refine_reftab[j][k - 1]);
				  }
				  nel.SetIndex(ind);
				  nel.flags.reverse = Refine_reverse[j];
				  if (elrev != 0)
				  {
				nel.flags.reverse = !nel.flags.reverse;
				swap(nel.PNum(3), nel.PNum(4));
				  }
        
				  if (j == 0)
				  {
		//C++ TO C# CONVERTER TODO TASK: The following line was determined to contain a copy constructor call - this should be verified and a copy constructor should be created:
		//ORIGINAL LINE: mesh.VolumeElement(ei) = nel;
					mesh.VolumeElement(new netgen.ElementIndex(ei)) = nel;
				  }
				  else
				  {
					mesh.AddVolumeElement(nel);
				  }
			   }
				break;
			  }
				  case ELEMENT_TYPE.HEX:
				  {
				 ArrayMem<PointIndex,27> pnums = new ArrayMem<PointIndex,27>(27);
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//		 static int betw[13][3] = { { 1, 2, 9 }, { 3, 4, 10 }, { 4, 1, 11 }, { 2, 3, 12 }, { 5, 6, 13 }, { 7, 8, 14 }, { 8, 5, 15 }, { 6, 7, 16 }, { 1, 5, 17 }, { 2, 6, 18 }, { 3, 7, 19 }, { 4, 8, 20 }, { 2, 8, 21 }};
        
					 /*
				 static int fbetw[12][3] =
				 { { 1, 3, 22 },
				   { 2, 4, 22 },
				   { 5, 7, 23 },
					   { 6, 8, 23 },
				   { 1, 6, 24 },
				   { 2, 5, 24 },
				   { 2, 7, 25 },
				   { 3, 6, 25 },
				   { 3, 8, 26 },
				   { 4, 7, 26 },
				   { 1, 8, 27 },
				   { 4, 5, 27 },
				   };
					 */
        
					 // updated by anonymous supporter, donations please to Karo W.
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//			 static int fbetw[12][3] = { { 11, 12, 22 }, { 9, 10, 22 }, { 13, 14, 23 }, { 15, 16, 23 }, { 9, 13, 24 }, { 17, 18, 24 }, { 12, 16, 25 }, { 18, 19, 25 }, { 19, 20, 26 }, { 10, 14, 26 }, { 11, 15, 27 }, { 17, 20, 27 }};
        
				 pnums = new PointIndex(-1);
        
				 for (int j = 1; j <= 8; j++)
				 {
					   pnums.Elem(j) = el.PNum(j);
				 }
        
        
				 for (int j = 0; j < 13; j++)
				 {
				   INDEX_2 i2 = new INDEX_2();
				   i2.I1() = pnums.Get(Refine_betw[j][0]);
				   i2.I2() = pnums.Get(Refine_betw[j][1]);
				   i2.Sort();
        
				   if (between.Used(i2))
				   {
					  pnums.Elem(9 + j) = between.Get(i2);
				   }
				   else
				   {
				  pnums.Elem(9 + j) = mesh.AddPoint(Center(new mesh.Point(i2.I1()), new mesh.Point(i2.I2())));
				  between.Set(i2, pnums.Elem(9 + j));
				   }
				 }
        
				for (int j = 0; j < 6; j++)
				{
				   INDEX_2 i2a = new INDEX_2();
				   INDEX_2 i2b = new INDEX_2();
				   i2a.I1() = pnums.Get(Refine_fbetw[2 * j][0]);
				   i2a.I2() = pnums.Get(Refine_fbetw[2 * j][1]);
				   i2a.Sort();
				   i2b.I1() = pnums.Get(Refine_fbetw[2 * j + 1][0]);
				   i2b.I2() = pnums.Get(Refine_fbetw[2 * j + 1][1]);
				   i2b.Sort();
        
				   if (between.Used(i2a))
				   {
				 pnums.Elem(22 + j) = between.Get(i2a);
				   }
				   else if (between.Used(i2b))
				   {
				 pnums.Elem(22 + j) = between.Get(i2b);
				   }
				   else
				   {
				   pnums.Elem(22 + j) = mesh.AddPoint(Center(new mesh.Point(i2a.I1()), new mesh.Point(i2a.I2())));
        
				   between.Set(i2a, pnums.Elem(22 + j));
				   }
				}
        
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//		static int reftab[8][8] = { { 1, 9, 22, 11, 17, 24, 21, 27 }, { 9, 2, 12, 22, 24, 18, 25, 21 }, { 11, 22, 10, 4, 27, 21, 26, 20}, { 22, 12, 3, 10, 21, 25, 19, 26}, { 17, 24, 21, 27, 5, 13, 23, 15}, { 24, 18, 25, 21, 13, 6, 16, 23}, { 27, 21, 26, 20, 15, 23, 14, 8}, { 21, 25, 19, 26, 23, 16, 7, 14} };
        
        
			   int ind = el.GetIndex();
			   for (int j = 0; j < 8; j++)
			   {
				  Element nel = new Element(ELEMENT_TYPE.HEX);
				  for (int k = 1; k <= 8; k++)
				  {
					nel.PNum(k) = pnums.Get(Refine_reftab[j][k - 1]);
				  }
				  nel.SetIndex(ind);
        
					  if (j == 0)
					  {
		//C++ TO C# CONVERTER TODO TASK: The following line was determined to contain a copy constructor call - this should be verified and a copy constructor should be created:
		//ORIGINAL LINE: mesh.VolumeElement(ei) = nel;
					mesh.VolumeElement(new netgen.ElementIndex(ei)) = nel;
					  }
				  else
				  {
					mesh.AddVolumeElement(nel);
				  }
			   }
				   break;
				  }
			  case ELEMENT_TYPE.PRISM:
			  {
				 ArrayMem<PointIndex,18> pnums = new ArrayMem<PointIndex,18>(18);
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//		 static int betw[9][3] = { { 3, 1, 7 }, { 1, 2, 8 }, { 3, 2, 9 }, { 6, 4, 10 }, { 4, 5, 11 }, { 6, 5, 12 }, { 1, 4, 13 }, { 3, 6, 14 }, { 2, 5, 15 }};
        
		// he: 15.jul 08, old version is wrong
		//                produces double points ad quad faces and inconsistent mesh
		// 	     static int fbetw[6][3] =
		// 	     { { 1, 6, 16 },
		// 	       { 3, 4, 16 },
		// 	       { 1, 5, 17 },
		//                { 2, 4, 17 },
		// 	       { 2, 6, 18 },
		// 	       { 3, 5, 18 },
		// 	       };
        
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//		   static int fbetw[6][3] = { { 7, 10, 16 }, { 14, 13, 16 }, { 11, 8, 17 }, { 13, 15, 17 }, { 12, 9, 18 }, { 14, 15, 18 }};
        
				 //int elrev = el.flags.reverse;
				   pnums = new PointIndex(-1);
        
				   for (int j = 1; j <= 6; j++)
				   {
				 pnums.Elem(j) = el.PNum(j);
				   }
				// if (elrev)
				// swap (pnums.Elem(3), pnums.Elem(4));
        
				 for (int j = 0; j < 9; j++)
				 {
				   INDEX_2 i2 = new INDEX_2();
				   i2.I1() = pnums.Get(Refine_betw[j][0]);
				   i2.I2() = pnums.Get(Refine_betw[j][1]);
				   i2.Sort();
        
				   if (between.Used(i2))
				   {
					  pnums.Elem(7 + j) = between.Get(i2);
				   }
				   else
				   {
				  pnums.Elem(7 + j) = mesh.AddPoint(Center(new mesh.Point(i2.I1()), new mesh.Point(i2.I2())));
				  between.Set(i2, pnums.Elem(7 + j));
				   }
				 }
        
				for (int j = 0; j < 3; j++)
				{
				   INDEX_2 i2a = new INDEX_2();
				   INDEX_2 i2b = new INDEX_2();
				   i2a.I1() = pnums.Get(Refine_fbetw[2 * j][0]);
				   i2a.I2() = pnums.Get(Refine_fbetw[2 * j][1]);
				   i2a.Sort();
				   i2b.I1() = pnums.Get(Refine_fbetw[2 * j + 1][0]);
				   i2b.I2() = pnums.Get(Refine_fbetw[2 * j + 1][1]);
				   i2b.Sort();
        
				   if (between.Used(i2a))
				   {
				 pnums.Elem(16 + j) = between.Get(i2a);
				   }
				   else if (between.Used(i2b))
				   {
				 pnums.Elem(16 + j) = between.Get(i2b);
				   }
				   else
				   {
				   pnums.Elem(16 + j) = mesh.AddPoint(Center(new mesh.Point(i2a.I1()), new mesh.Point(i2a.I2())));
        
				   between.Set(i2a, pnums.Elem(16 + j));
				   }
				}
        
        
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//		static int reftab[8][6] = { { 1, 8, 7, 13, 17, 16 }, { 7, 8, 9, 16, 17, 18 }, { 7, 9, 3, 16, 18, 14 }, { 8, 2, 9, 17, 15, 18 }, { 13, 17, 16, 4, 11, 10 }, { 16, 17, 18, 10, 11, 12 }, { 16, 18, 14, 10, 12, 6 }, { 17, 15, 18, 11, 5, 12 } };
        
        
			   int ind = el.GetIndex();
			   for (int j = 0; j < 8; j++)
			   {
				  Element nel = new Element(ELEMENT_TYPE.PRISM);
				  for (int k = 1; k <= 6; k++)
				  {
					nel.PNum(k) = pnums.Get(Refine_reftab[j][k - 1]);
				  }
				  nel.SetIndex(ind);
        
        
				  //nel.flags.reverse = reverse[j];
				  //if (elrev)
				 // {
				//nel.flags.reverse = 1 - nel.flags.reverse;
				//swap (nel.PNum(3), nel.PNum(4));
        
        
				  if (j == 0)
				  {
		//C++ TO C# CONVERTER TODO TASK: The following line was determined to contain a copy constructor call - this should be verified and a copy constructor should be created:
		//ORIGINAL LINE: mesh.VolumeElement(ei) = nel;
					mesh.VolumeElement(new netgen.ElementIndex(ei)) = nel;
				  }
				  else
				  {
					mesh.AddVolumeElement(nel);
				  }
			   }
				   break;
			  }
			  default:
				PrintSysError("Refine: undefined volume element type ", (int)el.GetType());
				break;
			}
			}
        
        
			// update identification tables
			for (int i = 1; i <= mesh.GetIdentifications().GetMaxNr(); i++)
			{
			Array<int,PointIndex.BASE> identmap = new Array<int,PointIndex.BASE>();
			mesh.GetIdentifications().GetMap(i, identmap);
        
			for (int j = 1; j <= between.GetNBags(); j++)
			{
			  for (int k = 1; k <= between.GetBagSize(j); k++)
			  {
				  INDEX_2 i2 = new INDEX_2();
				  PointIndex newpi = new PointIndex();
				  between.GetData(j, k, ref i2, ref newpi);
				  INDEX_2 oi2 = new INDEX_2(identmap.Get(i2.I1()), identmap.Get(i2.I2()));
				  oi2.Sort();
				  if (between.Used(oi2))
				  {
				  PointIndex onewpi = between.Get(oi2);
				  mesh.GetIdentifications().Add(new netgen.PointIndex(newpi), new netgen.PointIndex(onewpi), i);
				  }
			  }
			}
        
			}
        
			PrintMessage(5, "have 3d elements");
			mesh.ComputeNVertices();
			mesh.RebuildSurfaceElementLists();
			PrintMessage(5, "mesh updates complete");
			return;
        
			int cnttrials = 10;
			int wrongels = 0;
			for (int i = 1; i <= mesh.GetNE(); i++)
			{
			  if (mesh.VolumeElement(i).Volume(mesh.Points()) < 0)
			  {
			  wrongels++;
			  mesh.VolumeElement(i).flags.badel = 1;
			  }
			  else
			  {
			mesh.VolumeElement(i).flags.badel = 0;
			  }
			}
        
			if (wrongels != 0)
			{
			Console.Write("WARNING: ");
			Console.Write(wrongels);
			Console.Write(" with wrong orientation found");
			Console.Write("\n");
        
			int np = mesh.GetNP();
			Array<Point < 3>> should = new Array<Point < 3>>(np);
			Array<Point < 3>> can = new Array<Point < 3>>(np);
			for (int i = 1; i <= np; i++)
			{
				should.Elem(i) = can.Elem(i) = new mesh.Point(i);
			}
			for (int i = 1; i <= between.GetNBags(); i++)
			{
			  for (int j = 1; j <= between.GetBagSize(i); j++)
			  {
				  INDEX_2 parent = new INDEX_2();
				  PointIndex child = new PointIndex();
				  between.GetData(i, j, ref parent, ref child);
				  can.Elem(child) = Center(can.Elem(parent.I1()), can.Elem(parent.I2()));
			  }
			}
        
			BitArray boundp = new BitArray(np);
			boundp.Clear();
			foreach (var sel in mesh.SurfaceElements())
			{
				  foreach (var pi in sel.PNums())
				  {
					boundp.Set(pi);
				  }
			}
        
        
			double lam = 0.5;
        
			while (lam < 0.9 && cnttrials > 0)
			{
				lam = 2;
				do
				{
				lam *= 0.5;
				cnttrials--;
        
				Console.Write("lam = ");
				Console.Write(lam);
				Console.Write("\n");
        
				for (int i = 1; i <= np; i++)
				{
				  if (boundp.Test(i))
				  {
					  for (int j = 0; j < 3; j++)
					  {
					mesh.Point(i)(j) = lam * should.Get(i)(j) + (1 - lam) * can.Get(i)(j);
					  }
				  }
				  else
				  {
					mesh.Point(i) = can.Get(i);
				  }
				}
        
        
				BitArray free = new BitArray(mesh.GetNP());
				BitArray fhelp = new BitArray(mesh.GetNP());
				free.Clear();
				for (int i = 1; i <= mesh.GetNE(); i++)
				{
					Element el = mesh.VolumeElement(i);
					if (el.Volume(mesh.Points()) < 0)
					{
					  for (int j = 1; j <= el.GetNP(); j++)
					  {
					free.Set(el.PNum(j));
					  }
					}
				}
				for (int k = 1; k <= 3; k++)
				{
					fhelp.Clear();
					for (int i = 1; i <= mesh.GetNE(); i++)
					{
					Element el = mesh.VolumeElement(i);
					int freeel = 0;
					for (int j = 1; j <= el.GetNP(); j++)
					{
					  if (free.Test(el.PNum(j)))
					  {
						freeel = 1;
					  }
					}
					if (freeel != 0)
					{
					  for (int j = 1; j <= el.GetNP(); j++)
					  {
						fhelp.Set(el.PNum(j));
					  }
					}
					}
					free.Or(fhelp);
				}
        
				(*testout) << "smooth points: " << "\n";
				for (int i = 1; i <= free.Size(); i++)
				{
				  if (free.Test(i))
				  {
					(*testout) << "p " << i << "\n";
				  }
				}
        
				(*testout) << "surf points: " << "\n";
				foreach (var sel in mesh.SurfaceElements())
				{
				  foreach (var pi in sel.PNums())
				  {
					(*testout) << pi << "\n";
				  }
				}
        
				mesh.CalcSurfacesOfNode();
				free.Invert();
				mesh.FixPoints(free);
				MeshingParameters dummymp = new MeshingParameters();
				mesh.ImproveMesh(dummymp, OPTIMIZEGOAL.OPT_REST);
        
        
				wrongels = 0;
				for (int i = 1; i <= mesh.GetNE(); i++)
				{
					if (mesh.VolumeElement(i).Volume(mesh.Points()) < 0)
					{
					wrongels++;
					mesh.VolumeElement(i).flags.badel = 1;
					(*testout) << "wrong el: ";
					for (int j = 1; j <= 4; j++)
					{
					  (*testout) << mesh.VolumeElement(i).PNum(j) << " ";
					}
					(*testout) << "\n";
					}
					else
					{
					  mesh.VolumeElement(i).flags.badel = 0;
					}
				}
				Console.Write("wrongels = ");
				Console.Write(wrongels);
				Console.Write("\n");
				} while (wrongels != 0 && cnttrials > 0);
        
				for (int i = 1; i <= np; i++)
				{
				  can.Elem(i) = new mesh.Point(i);
				}
			}
			}
        
			if (cnttrials <= 0)
			{
			cerr << "ERROR: Sorry, reverted elements" << "\n";
			}
        
			mesh.ComputeNVertices();
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void MakeSecondOrder(Mesh mesh)
		  {
		//C++ TO C# CONVERTER TODO TASK: There is no equivalent to 'const_cast' in C#:
			const_cast<Refinement&> (this).MakeSecondOrder(mesh);
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void MakeSecondOrder(Mesh mesh)
		  {
			/*
			  Berlin, 2014: if we have curved surface elements, keep them !
			*/
        
			mesh.ComputeNVertices();
			// mesh.SetNP(mesh.GetNV());
			mesh.SetNP(mesh.GetNP()); // setup multilevel-table
        
			INDEX_2_HASHTABLE<PointIndex> between = new INDEX_2_HASHTABLE<PointIndex>(mesh.GetNP() + 5);
        
			for (SurfaceElementIndex sei = 0; sei < mesh.GetNSE(); sei++)
			{
			Element2d el = mesh[sei];
        
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static int betw_trig[3][3] = { { 1, 2, 3 }, { 0, 2, 4 }, { 0, 1, 5 } };
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static int betw_quad6[2][3] = { { 0, 1, 4 }, { 3, 2, 5 } };
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static int betw_quad8[4][3] = { { 0, 1, 4 }, { 3, 2, 5 }, { 0, 3, 6 }, { 1, 2, 7 } };
        
			int onp = 0;
			int[] betw = 0;
			switch (el.GetType())
			{
			  case ELEMENT_TYPE.TRIG6:
			  {
				  betw = MakeSecondOrder_betw_trig;
				  onp = 3;
				  break;
			  }
			  case ELEMENT_TYPE.QUAD6:
			  {
					  betw = MakeSecondOrder_betw_quad6;
					  onp = 4;
					  break;
			  }
			  case ELEMENT_TYPE.QUAD8:
			  {
					  betw = MakeSecondOrder_betw_quad8;
					  onp = 4;
					  break;
			  }
				  default:
					;
					break;
			}
        
				if (betw)
				{
				  for (int j = 0; j < el.GetNP() - onp; j++)
				  {
					  int pi1 = el[betw[j][0]];
					  int pi2 = el[betw[j][1]];
					  INDEX_2 i2 = INDEX_2.Sort(pi1, pi2);
					  between.Set(i2, el[onp + j]);
				  }
				}
			}
        
        
			bool thinlayers = false;
			for (ElementIndex ei = 0; ei < mesh.GetNE(); ei++)
			{
			  if (mesh[ei].GetType() == ELEMENT_TYPE.PRISM || mesh[ei].GetType() == ELEMENT_TYPE.PRISM12)
			  {
			thinlayers = true;
			  }
			}
        
        
			int nseg = mesh.GetNSeg();
			for (SegmentIndex si = 0; si < nseg; si++)
			{
		//C++ TO C# CONVERTER TODO TASK: The following line was determined to contain a copy constructor call - this should be verified and a copy constructor should be created:
		//ORIGINAL LINE: Segment & el = mesh.LineSegment(si);
			Segment el = mesh.LineSegment(new netgen.SegmentIndex(si));
        
			INDEX_2 i2 = INDEX_2.Sort(new netgen.Segment(el[0]), new netgen.Segment(el[1]));
        
			if (between.Used(i2))
			{
			  el[2] = between.Get(i2);
			}
			else
			{
				Point < 3> pb;
				EdgePointGeomInfo ngi = new EdgePointGeomInfo();
					PointBetween(new mesh.Point(el[0]), new mesh.Point(el[1]), 0.5, el.surfnr1, el.surfnr2, el.epgeominfo[0], el.epgeominfo[1], pb, ngi);
        
				el[2] = mesh.AddPoint(pb, new mesh.Point(el[0]).GetLayer(), POINTTYPE.EDGEPOINT);
				between.Set(i2, el[2]);
			}
			}
        
			// refine surface elements
			for (SurfaceElementIndex sei = 0; sei < mesh.GetNSE(); sei++)
			{
			int j;
		//C++ TO C# CONVERTER TODO TASK: The following line was determined to contain a copy constructor call - this should be verified and a copy constructor should be created:
		//ORIGINAL LINE: const Element2d & el = mesh.SurfaceElement(sei);
			Element2d el = mesh.SurfaceElement(new netgen.SurfaceElementIndex(sei));
        
			int onp = 0;
        
			Element2d newel = new Element2d(ELEMENT_TYPE.TRIG);
			newel.SetIndex(el.GetIndex());
        
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static int betw_trig[3][3] = { { 1, 2, 3 }, { 0, 2, 4 }, { 0, 1, 5 } };
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static int betw_quad6[2][3] = { { 0, 1, 4 }, { 3, 2, 5 } };
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static int betw_quad8[4][3] = { { 0, 1, 4 }, { 3, 2, 5 }, { 0, 3, 6 }, { 1, 2, 7 } };
			int[] betw = 0;
        
			switch (el.GetType())
			{
			  case ELEMENT_TYPE.TRIG:
			  case ELEMENT_TYPE.TRIG6:
			  {
				  betw = MakeSecondOrder_betw_trig;
				  newel.SetType(ELEMENT_TYPE.TRIG6);
				  onp = 3;
				  break;
			  }
			  case ELEMENT_TYPE.QUAD:
			  case ELEMENT_TYPE.QUAD6:
			  case ELEMENT_TYPE.QUAD8:
			  {
				  if (thinlayers)
				  {
				  betw = MakeSecondOrder_betw_quad6;
				  newel.SetType(ELEMENT_TYPE.QUAD6);
				  }
				  else
				  {
				  betw = MakeSecondOrder_betw_quad8;
				  newel.SetType(ELEMENT_TYPE.QUAD8);
				  }
				  onp = 4;
				  break;
			  }
			  default:
				PrintSysError("Unhandled element in secondorder:", (int)el.GetType());
				break;
			}
        
			for (j = 0; j < onp; j++)
			{
			  newel[j] = el[j];
			}
        
			int nnp = newel.GetNP();
			for (j = 0; j < nnp - onp; j++)
			{
				int pi1 = newel[betw[j][0]];
				int pi2 = newel[betw[j][1]];
        
				INDEX_2 i2 = INDEX_2.Sort(pi1, pi2);
        
				if (between.Used(i2))
				{
				  newel[onp + j] = between.Get(i2);
				}
				else
				{
				Point < 3> pb;
				PointGeomInfo newgi = new PointGeomInfo();
				PointBetween(new mesh.Point(pi1), new mesh.Point(pi2), 0.5, mesh.GetFaceDescriptor(el.GetIndex()).SurfNr(), el.GeomInfoPi(betw[j][0] + 1), el.GeomInfoPi(betw[j][1] + 1), pb, newgi);
        
				newel[onp + j] = mesh.AddPoint(pb, new mesh.Point(pi1).GetLayer(), POINTTYPE.SURFACEPOINT);
				between.Set(i2, newel[onp + j]);
				}
			}
        
		//C++ TO C# CONVERTER TODO TASK: The following line was determined to contain a copy constructor call - this should be verified and a copy constructor should be created:
		//ORIGINAL LINE: mesh.SurfaceElement(sei) = newel;
			mesh.SurfaceElement(new netgen.SurfaceElementIndex(sei)) = newel;
			}
        
        
			//    int i, j;
        
        
        
			// refine volume elements
			for (int i = 1; i <= mesh.GetNE(); i++)
			{
			Element el = mesh.VolumeElement(i);
			int onp = 0;
        
			Element newel = new Element(ELEMENT_TYPE.TET);
			newel.SetIndex(el.GetIndex());
        
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static int betw_tet[6][3] = { { 0, 1, 4 }, { 0, 2, 5 }, { 0, 3, 6 }, { 1, 2, 7 }, { 1, 3, 8 }, { 2, 3, 9 } };
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static int betw_prism[6][3] = { { 0, 2, 6 }, { 0, 1, 7 }, { 1, 2, 8 }, { 3, 5, 9 }, { 3, 4, 10 }, { 4, 5, 11 }};
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//		static int betw_prism15[9][3] = { { 0, 1, 6 }, { 0, 2, 7 }, { 1, 2, 8 }, { 0, 3, 9 }, { 1, 4, 10 }, { 2, 5, 11 }, { 3, 4, 12 }, { 3, 5, 13 }, { 4, 5, 14 } };
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//		static int betw_pyramid[8][3] = { { 0, 1, 5 }, { 3, 2, 6 }, { 3, 0, 7 }, { 1, 2, 8 }, { 0, 4, 9 }, { 1, 4, 10 }, { 2, 4, 11 }, { 3, 4, 12 } };
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static int betw_hex[12][3] = { { 0, 1, 8 }, { 2, 3, 9 }, { 3, 0, 10 }, { 1, 2, 11 }, { 4, 5, 12 }, { 6, 7, 13 }, { 7, 4, 14 }, { 5, 6, 15 }, { 0, 4, 16 }, { 1, 5, 17 }, { 2, 6, 18 }, { 3, 7, 19 }};
        
			int[] betw = 0;
        
			switch (el.GetType())
			{
			  case ELEMENT_TYPE.TET:
			  case ELEMENT_TYPE.TET10:
			  {
				  betw = MakeSecondOrder_betw_tet;
				  newel.SetType(ELEMENT_TYPE.TET10);
				  onp = 4;
				  break;
			  }
			  case ELEMENT_TYPE.PRISM:
			  case ELEMENT_TYPE.PRISM12:
			  {
				  betw = MakeSecondOrder_betw_prism;
				  newel.SetType(ELEMENT_TYPE.PRISM12);
				  onp = 6;
				  break;
			  }
				  case ELEMENT_TYPE.PRISM15:
				  {
					  betw = MakeSecondOrder_betw_prism15;
					  newel.SetType(ELEMENT_TYPE.PRISM15);
					  onp = 6;
					  break;
				  }
				  case ELEMENT_TYPE.PYRAMID:
				  case ELEMENT_TYPE.PYRAMID13:
				  {
					  betw = MakeSecondOrder_betw_pyramid;
					  newel.SetType(ELEMENT_TYPE.PYRAMID13);
					  onp = 5;
					  break;
				  }
				  case ELEMENT_TYPE.HEX:
				  case ELEMENT_TYPE.HEX20:
				  {
				  betw = MakeSecondOrder_betw_hex;
				  newel.SetType(ELEMENT_TYPE.HEX20);
				  onp = 8;
				  break;
				  }
			  default:
				PrintSysError("MakeSecondOrder, illegal vol type ", (int)el.GetType());
				break;
			}
        
        
			for (int j = 1; j <= onp; j++)
			{
			  newel.PNum(j) = el.PNum(j);
			}
			int nnp = newel.GetNP();
        
			for (int j = 0; j < nnp - onp; j++)
			{
				INDEX_2 i2 = new INDEX_2(newel[betw[j][0]], newel[betw[j][1]]);
				i2.Sort();
        
				if (between.Used(i2))
				{
				  newel.PNum(onp + 1 + j) = between.Get(i2);
				}
				else
				{
				newel.PNum(onp + 1 + j) = mesh.AddPoint(Center(new mesh.Point(i2.I1()), new mesh.Point(i2.I2())), new mesh.Point(i2.I1()).GetLayer(), POINTTYPE.INNERPOINT);
        
				between.Set(i2, newel.PNum(onp + 1 + j));
				}
			}
        
			mesh.VolumeElement(i) = newel;
			}
        
        
			// makes problems after linear mesh refinement, since
			// 2nd order identifications are not removed
			// update identification tables
			for (int i = 1; i <= mesh.GetIdentifications().GetMaxNr(); i++)
			{
			Array<int,PointIndex.BASE> identmap = new Array<int,PointIndex.BASE>();
			mesh.GetIdentifications().GetMap(i, identmap);
        
			for (INDEX_2_HASHTABLE<PointIndex>.Iterator it = between.Begin(); it != between.End(); it++)
			{
				  INDEX_2 i2 = new INDEX_2();
				  PointIndex newpi = new PointIndex();
				  between.GetData(it, ref i2, ref newpi);
				  INDEX_2 oi2 = new INDEX_2(identmap.Get(i2.I1()), identmap.Get(i2.I2()));
				  oi2.Sort();
				  if (between.Used(oi2))
				  {
				  PointIndex onewpi = between.Get(oi2);
				  mesh.GetIdentifications().Add(new netgen.PointIndex(newpi), new netgen.PointIndex(onewpi), i);
				  }
			}
        
			/*
			for (int j = 1; j <= between.GetNBags(); j++)
			  for (int k = 1; k <= between.GetBagSize(j); k++)
			    {
			      INDEX_2 i2;
			      int newpi;
			      between.GetData (j, k, i2, newpi);
			      INDEX_2 oi2(identmap.Get(i2.I1()),
					  identmap.Get(i2.I2()));
			      oi2.Sort();
			      if (between.Used (oi2))
				{
				  int onewpi = between.Get(oi2);
				  mesh.GetIdentifications().Add (newpi, onewpi, i);
				}
			    }
			*/
			}
        
        
			//  mesh.mglevels++;
			int oldsize = mesh.mlbetweennodes.Size();
			mesh.mlbetweennodes.SetSize(mesh.GetNP());
			for (int i = oldsize; i < mesh.GetNP(); i++)
			{
			  mesh.mlbetweennodes[i] = new INDEX_2(0, 0);
			}
        
			/*
			for (i = 1; i <= between.GetNBags(); i++)
			  for (j = 1; j <= between.GetBagSize(i); j++)
			{
			  INDEX_2 oldp;
			  int newp;
			  between.GetData (i, j, oldp, newp);
			  mesh.mlbetweennodes.Elem(newp) = oldp;
			}
			*/
        
			for (INDEX_2_HASHTABLE<PointIndex>.Iterator it = between.Begin(); it != between.End(); it++)
			{
			mesh.mlbetweennodes[between.GetData(it)] = between.GetHash(it);
			}
        
			mesh.ComputeNVertices();
			mesh.RebuildSurfaceElementLists();
			//  ValidateSecondOrder (mesh);
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void ValidateSecondOrder(Mesh mesh)
		  {
			PrintMessage(3, "Validate mesh");
			int np = mesh.GetNP();
			int ne = mesh.GetNE();
			// int i, j;
			Array<INDEX_2> parents = new Array<INDEX_2>(np);
        
			for (int i = 1; i <= np; i++)
			{
			  parents.Elem(i) = new INDEX_2(0, 0);
			}
        
			for (int i = 1; i <= ne; i++)
			{
			Element el = mesh.VolumeElement(i);
			if (el.GetType() == ELEMENT_TYPE.TET10)
			{
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//		static int betweentab[6][3] = { { 1, 2, 5 }, { 1, 3, 6 }, { 1, 4, 7 }, { 2, 3, 8 }, { 2, 4, 9 }, { 3, 4, 10 } };
				for (int j = 0; j < 6; j++)
				{
				int f1 = el.PNum(ValidateSecondOrder_betweentab[j][0]);
				int f2 = el.PNum(ValidateSecondOrder_betweentab[j][1]);
				int son = el.PNum(ValidateSecondOrder_betweentab[j][2]);
				parents.Elem(son).I1() = f1;
				parents.Elem(son).I2() = f2;
				}
			}
			}
        
			ValidateRefinedMesh(mesh, parents);
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void ValidateRefinedMesh(Mesh mesh, Array<INDEX_2> parents)
		  {
			// int i, j, k;
        
			// homotopy method
        
			int ne = mesh.GetNE();
        
			int cnttrials = 100;
			int wrongels = 0;
			for (int i = 1; i <= ne; i++)
			{
			  if (mesh.VolumeElement(i).CalcJacobianBadness(mesh.Points()) > 1e10)
			  {
			  wrongels++;
			  mesh.VolumeElement(i).flags.badel = 1;
			  }
			  else
			  {
			mesh.VolumeElement(i).flags.badel = 0;
			  }
			}
        
			double facok = 0;
			double factry;
        
			BitArray illegalels = new BitArray(ne);
			illegalels.Clear();
        
        
			if (wrongels != 0)
			{
			Console.Write("WARNING: ");
			Console.Write(wrongels);
			Console.Write(" illegal element(s) found");
			Console.Write("\n");
        
			int np = mesh.GetNP();
			Array<Point < 3>> should = new Array<Point < 3>>(np);
			Array<Point < 3>> can = new Array<Point < 3>>(np);
        
			for (int i = 1; i <= np; i++)
			{
				should.Elem(i) = can.Elem(i) = new mesh.Point(i);
			}
        
			for (int i = 1; i <= parents.Size(); i++)
			{
				if (parents.Get(i).I1())
				{
				  can.Elem(i) = Center(can.Elem(parents.Get(i).I1()), can.Elem(parents.Get(i).I2()));
				}
			}
        
			BitArray boundp = new BitArray(np);
			boundp.Clear();
			for (int i = 1; i <= mesh.GetNSE(); i++)
			{
				Element2d sel = mesh.SurfaceElement(i);
				for (int j = 1; j <= sel.GetNP(); j++)
				{
				  boundp.Set(sel.PNum(j));
				}
			}
        
        
			(*testout) << "bpoints:" << "\n";
			for (int i = 1; i <= np; i++)
			{
			  if (boundp.Test(i))
			  {
				(*testout) << i << "\n";
			  }
			}
        
			double lam = 0.5;
        
			while (facok < 1 - 1e-8 && cnttrials > 0)
			{
				lam *= 4;
				if (lam > 2)
				{
					lam = 2;
				}
        
				do
				{
				//	      cout << "trials: " << cnttrials << endl;
				lam *= 0.5;
				cnttrials--;
        
				Console.Write("lam = ");
				Console.Write(lam);
				Console.Write("\n");
        
				factry = lam + (1 - lam) * facok;
				Console.Write("trying: ");
				Console.Write(factry);
				Console.Write("\n");
        
				for (int i = 1; i <= np; i++)
				{
				  if (boundp.Test(i))
				  {
					  for (int j = 0; j < 3; j++)
					  {
					mesh.Point(i)(j) = lam * should.Get(i)(j) + (1 - lam) * can.Get(i)(j);
					  }
				  }
				  else
				  {
					mesh.Point(i) = Point < 3> (can.Get(i));
				  }
				}
        
				//	      (*testout) << "bad els: " << endl;
				wrongels = 0;
				for (int i = 1; i <= ne; i++)
				{
					if (!illegalels.Test(i) && mesh.VolumeElement(i).CalcJacobianBadness(mesh.Points()) > 1e10)
					{
					wrongels++;
					Element el = mesh.VolumeElement(i);
					el.flags.badel = 1;
        
        
					if (lam < 1e-4)
					{
					  illegalels.Set(i);
					}
        
        
					/*
					  (*testout) << i << ": ";
					  for (j = 1; j <= el.GetNP(); j++)
					  (*testout) << el.PNum(j) << " ";
					  (*testout) << endl;
					*/
					}
					else
					{
					  mesh.VolumeElement(i).flags.badel = 0;
					}
				}
				Console.Write("wrongels = ");
				Console.Write(wrongels);
				Console.Write("\n");
				} while (wrongels != 0 && cnttrials > 0);
        
				mesh.CalcSurfacesOfNode();
				MeshingParameters dummymp = new MeshingParameters();
				mesh.ImproveMeshJacobian(dummymp, OPTIMIZEGOAL.OPT_WORSTCASE);
        
				facok = factry;
				for (int i = 1; i <= np; i++)
				{
				  can.Elem(i) = new mesh.Point(i);
				}
			}
			}
        
        
        
			for (int i = 1; i <= ne; i++)
			{
			if (illegalels.Test(i))
			{
				Console.Write("illegal element: ");
				Console.Write(i);
				Console.Write("\n");
				mesh.VolumeElement(i).flags.badel = 1;
			}
			else
			{
			  mesh.VolumeElement(i).flags.badel = 0;
			}
			}
        
			/*
			  if (cnttrials <= 0)
			  {
			  cerr << "ERROR: Sorry, illegal elements:" << endl;
			  }
			*/
		  }
	}
}