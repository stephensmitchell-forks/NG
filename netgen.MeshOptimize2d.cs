namespace netgen
{

	public class MeshOptimize2d
	{
//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void GenericImprove(Mesh mesh)
		  {
			if (!faceindex)
			{
			if (writestatus)
			{
			  PrintMessage(3, "Generic Improve");
			}
        
			for (faceindex = 1; faceindex <= mesh.GetNFD(); faceindex++)
			{
			  GenericImprove(mesh);
			}
        
			faceindex = 0;
			}
        
			// int j, k, l, ri;
			int np = mesh.GetNP();
			int ne = mesh.GetNSE();
			//    SurfaceElementIndex sei;
        
        
		//     for (SurfaceElementIndex sei = 0; sei < ne; sei++)
		//       {
		// 	const Element2d & el = mesh[sei];
		// 	(*testout) << "element " << sei << ": " <<flush;
		// 	for(int j=0; j<el.GetNP(); j++)
		// 	  (*testout) << el[j] << " " << flush;
		// 	(*testout) << "IsDeleted() " << el.IsDeleted()<< endl;
		//       }
        
			bool ok;
			int olddef;
			int newdef;
        
			Array<ImprovementRule> rules = new Array<ImprovementRule>();
			Array<SurfaceElementIndex> elmap = new Array<SurfaceElementIndex>();
			Array<int> elrot = new Array<int>();
			Array<PointIndex> pmap = new Array<PointIndex>();
			Array<PointGeomInfo> pgi = new Array<PointGeomInfo>();
        
			int surfnr = mesh.GetFaceDescriptor(faceindex).SurfNr();
        
			ImprovementRule r1;
        
			// 2 triangles to quad
			r1 = new ImprovementRule();
			r1.oldels.Append(new Element2d(1, 2, 3));
			r1.oldels.Append(new Element2d(3, 2, 4));
			r1.newels.Append(new Element2d(1, 2, 4, 3));
			r1.deledges.Append(new INDEX_2(2, 3));
			r1.onp = 4;
			r1.bonus = 2;
			rules.Append(r1);
        
			// 2 quad to 1 quad
			r1 = new ImprovementRule();
			r1.oldels.Append(new Element2d(1, 2, 3, 4));
			r1.oldels.Append(new Element2d(4, 3, 2, 5));
			r1.newels.Append(new Element2d(1, 2, 5, 4));
			r1.deledges.Append(new INDEX_2(2, 3));
			r1.deledges.Append(new INDEX_2(3, 4));
			r1.onp = 5;
			r1.bonus = 0;
			rules.Append(r1);
        
			// swap quads
			r1 = new ImprovementRule();
			r1.oldels.Append(new Element2d(1, 2, 3, 4));
			r1.oldels.Append(new Element2d(3, 2, 5, 6));
			r1.newels.Append(new Element2d(1, 6, 3, 4));
			r1.newels.Append(new Element2d(1, 2, 5, 6));
			r1.deledges.Append(new INDEX_2(2, 3));
			r1.onp = 6;
			r1.bonus = 0;
			rules.Append(r1);
        
			// three quads to 2
			r1 = new ImprovementRule();
			r1.oldels.Append(new Element2d(1, 2, 3, 4));
			r1.oldels.Append(new Element2d(2, 5, 6, 3));
			r1.oldels.Append(new Element2d(3, 6, 7, 4));
			r1.newels.Append(new Element2d(1, 2, 5, 4));
			r1.newels.Append(new Element2d(4, 5, 6, 7));
			r1.deledges.Append(new INDEX_2(2, 3));
			r1.deledges.Append(new INDEX_2(3, 4));
			r1.deledges.Append(new INDEX_2(3, 6));
			r1.onp = 7;
			r1.bonus = -1;
			rules.Append(r1);
        
			// quad + 2 connected trigs to quad
			r1 = new ImprovementRule();
			r1.oldels.Append(new Element2d(1, 2, 3, 4));
			r1.oldels.Append(new Element2d(2, 5, 3));
			r1.oldels.Append(new Element2d(3, 5, 4));
			r1.newels.Append(new Element2d(1, 2, 5, 4));
			r1.deledges.Append(new INDEX_2(2, 3));
			r1.deledges.Append(new INDEX_2(3, 4));
			r1.deledges.Append(new INDEX_2(3, 5));
			r1.onp = 5;
			r1.bonus = 0;
			rules.Append(r1);
        
			// quad + 2 non-connected trigs to quad (a and b)
			r1 = new ImprovementRule();
			r1.oldels.Append(new Element2d(1, 2, 3, 4));
			r1.oldels.Append(new Element2d(2, 6, 3));
			r1.oldels.Append(new Element2d(1, 4, 5));
			r1.newels.Append(new Element2d(1, 3, 4, 5));
			r1.newels.Append(new Element2d(1, 2, 6, 3));
			r1.deledges.Append(new INDEX_2(1, 4));
			r1.deledges.Append(new INDEX_2(2, 3));
			r1.onp = 6;
			r1.bonus = 0;
			rules.Append(r1);
        
			r1 = new ImprovementRule();
			r1.oldels.Append(new Element2d(1, 2, 3, 4));
			r1.oldels.Append(new Element2d(2, 6, 3));
			r1.oldels.Append(new Element2d(1, 4, 5));
			r1.newels.Append(new Element2d(1, 2, 4, 5));
			r1.newels.Append(new Element2d(4, 2, 6, 3));
			r1.deledges.Append(new INDEX_2(1, 4));
			r1.deledges.Append(new INDEX_2(2, 3));
			r1.onp = 6;
			r1.bonus = 0;
			rules.Append(r1);
        
			// two quad + trig -> one quad + trig
			r1 = new ImprovementRule();
			r1.oldels.Append(new Element2d(1, 2, 3, 4));
			r1.oldels.Append(new Element2d(2, 5, 6, 3));
			r1.oldels.Append(new Element2d(4, 3, 6));
			r1.newels.Append(new Element2d(1, 2, 6, 4));
			r1.newels.Append(new Element2d(2, 5, 6));
			r1.deledges.Append(new INDEX_2(2, 3));
			r1.deledges.Append(new INDEX_2(3, 4));
			r1.deledges.Append(new INDEX_2(3, 6));
			r1.onp = 6;
			r1.bonus = -1;
			rules.Append(r1);
        
			// swap quad + trig (a and b)
			r1 = new ImprovementRule();
			r1.oldels.Append(new Element2d(1, 2, 3, 4));
			r1.oldels.Append(new Element2d(2, 5, 3));
			r1.newels.Append(new Element2d(2, 5, 3, 4));
			r1.newels.Append(new Element2d(1, 2, 4));
			r1.deledges.Append(new INDEX_2(2, 3));
			r1.onp = 5;
			r1.bonus = 0;
			rules.Append(r1);
        
			r1 = new ImprovementRule();
			r1.oldels.Append(new Element2d(1, 2, 3, 4));
			r1.oldels.Append(new Element2d(2, 5, 3));
			r1.newels.Append(new Element2d(1, 2, 5, 3));
			r1.newels.Append(new Element2d(1, 3, 4));
			r1.deledges.Append(new INDEX_2(2, 3));
			r1.onp = 5;
			r1.bonus = 0;
			rules.Append(r1);
        
        
			// 2 quads to quad + 2 trigs
			r1 = new ImprovementRule();
			r1.oldels.Append(new Element2d(1, 2, 3, 4));
			r1.oldels.Append(new Element2d(3, 2, 5, 6));
			r1.newels.Append(new Element2d(1, 5, 6, 4));
			r1.newels.Append(new Element2d(1, 2, 5));
			r1.newels.Append(new Element2d(4, 6, 3));
			r1.deledges.Append(new INDEX_2(2, 3));
			r1.onp = 6;
			r1.bonus = 0;
			//    rules.Append (r1);
        
        
        
        
			Array<int> mapped = new Array<int>(rules.Size());
			Array<int> used = new Array<int>(rules.Size());
			used = 0;
			mapped = 0;
        
        
        
			for (int ri = 0; ri < rules.Size(); ri++)
			{
			ImprovementRule rule = *rules[ri];
			rule.incelsonnode.SetSize(rule.onp);
			rule.reused.SetSize(rule.onp);
        
			for (int j = 1; j <= rule.onp; j++)
			{
				rule.incelsonnode.Elem(j) = 0;
				rule.reused.Elem(j) = 0;
			}
        
			for (int j = 1; j <= rule.oldels.Size(); j++)
			{
				Element2d el = rule.oldels.Elem(j);
				for (int k = 1; k <= el.GetNP(); k++)
				{
				  rule.incelsonnode.Elem(el.PNum(k))--;
				}
			}
        
			for (int j = 1; j <= rule.newels.Size(); j++)
			{
				Element2d el = rule.newels.Elem(j);
				for (int k = 1; k <= el.GetNP(); k++)
				{
				rule.incelsonnode.Elem(el.PNum(k))++;
				rule.reused.Elem(el.PNum(k)) = 1;
				}
			}
			}
        
        
        
        
			TABLE<int,PointIndex.BASE> elonnode = new TABLE<int,PointIndex.BASE>(np);
			Array<int,PointIndex.BASE> nelonnode = new Array<int,PointIndex.BASE>(np);
			TABLE<SurfaceElementIndex> nbels = new TABLE<SurfaceElementIndex>(ne);
        
			nelonnode = -4;
        
			for (SurfaceElementIndex sei = 0; sei < ne; sei++)
			{
			Element2d el = mesh[sei];
        
			if (el.GetIndex() == faceindex && !el.IsDeleted())
			{
				for (int j = 0; j < el.GetNP(); j++)
				{
				  elonnode.Add(el[j], sei);
				}
			}
			if (!el.IsDeleted())
			{
				for (int j = 0; j < el.GetNP(); j++)
				{
				  nelonnode[el[j]]++;
				}
			}
			}
        
			for (SurfaceElementIndex sei = 0; sei < ne; sei++)
			{
			Element2d el = mesh[sei];
			if (el.GetIndex() == faceindex && !el.IsDeleted())
			{
				for (int j = 0; j < el.GetNP(); j++)
				{
				for (int k = 0; k < elonnode[el[j]].Size(); k++)
				{
					int nbel = elonnode[el[j]][k];
					bool inuse = false;
					for (int l = 0; l < nbels[sei].Size(); l++)
					{
					  if (nbels[sei][l] == nbel)
					  {
					inuse = true;
					  }
					}
					if (!inuse)
					{
					  nbels.Add(sei, nbel);
					}
				}
				}
			}
			}
        
        
			for (int ri = 0; ri < rules.Size(); ri++)
			{
			ImprovementRule rule = *rules[ri];
        
			elmap.SetSize(rule.oldels.Size());
			elrot.SetSize(rule.oldels.Size());
			pmap.SetSize(rule.onp);
			pgi.SetSize(rule.onp);
        
        
			for (SurfaceElementIndex sei = 0; sei < ne; sei++)
			{
				if (multithread.terminate)
				{
				  break;
				}
				if (mesh[sei].IsDeleted())
				{
					continue;
				}
        
				elmap[0] = sei;
				FlatArray<SurfaceElementIndex> neighbours = nbels[sei];
        
				for (elrot[0] = 0; elrot[0] < mesh[sei].GetNP(); elrot[0]++)
				{
				Element2d el0 = mesh[sei];
				Element2d rel0 = rule.oldels[0];
        
				if (el0.GetIndex() != faceindex)
				{
					continue;
				}
				if (el0.IsDeleted())
				{
					continue;
				}
				if (el0.GetNP() != rel0.GetNP())
				{
					continue;
				}
        
        
				pmap = new PointIndex(-1);
        
				for (int k = 0; k < el0.GetNP(); k++)
				{
					pmap.Elem(rel0[k]) = el0.PNumMod(k + elrot[0] + 1);
					pgi.Elem(rel0[k]) = el0.GeomInfoPiMod(k + elrot[0] + 1);
				}
        
				ok = true;
				for (int i = 1; i < elmap.Size(); i++)
				{
					// try to find a mapping for reference-element i
        
					Element2d rel = rule.oldels[i];
					bool possible = false;
        
					for (elmap[i] = 0; elmap[i] < neighbours.Size(); elmap[i]++)
					{
					Element2d el = mesh[neighbours[elmap[i]]];
					if (el.IsDeleted())
					{
						continue;
					}
					if (el.GetNP() != rel.GetNP())
					{
						continue;
					}
        
					for (elrot[i] = 0; elrot[i] < rel.GetNP(); elrot[i]++)
					{
						possible = true;
        
						for (int k = 0; k < rel.GetNP(); k++)
						{
						  if (pmap.Elem(rel[k]) != -1 && pmap.Elem(rel[k]) != el.PNumMod(k + elrot[i] + 1))
						  {
						possible = false;
						  }
						}
        
						if (possible)
						{
						for (int k = 0; k < el.GetNP(); k++)
						{
							pmap.Elem(rel[k]) = el.PNumMod(k + elrot[i] + 1);
							pgi.Elem(rel[k]) = el.GeomInfoPiMod(k + elrot[i] + 1);
						}
						break;
						}
					}
					if (possible)
					{
						break;
					}
					}
        
					if (!possible)
					{
					ok = false;
					break;
					}
        
					elmap[i] = neighbours[elmap[i]];
				}
        
				for (int i = 0; ok && i < rule.deledges.Size(); i++)
				{
					ok = !mesh.IsSegment(pmap.Elem(rule.deledges[i].I1()), pmap.Elem(rule.deledges[i].I2()));
				}
        
        
        
        
				if (!ok)
				{
					continue;
				}
        
				mapped[ri]++;
        
				olddef = 0;
				for (int j = 1; j <= pmap.Size(); j++)
				{
				  olddef += sqr(nelonnode[pmap.Get(j)]);
				}
				olddef += rule.bonus;
        
				newdef = 0;
				for (int j = 1; j <= pmap.Size(); j++)
				{
				  if (rule.reused.Get(j))
				  {
					newdef += sqr(nelonnode[pmap.Get(j)] + rule.incelsonnode.Get(j));
				  }
				}
        
				if (newdef > olddef)
				{
				  continue;
				}
        
				// calc metric badness
				double bad1 = 0;
				double bad2 = 0;
				Vec < 3> n;
        
				SelectSurfaceOfPoint(new mesh.Point(pmap.Get(1)), pgi.Get(1));
				GetNormalVector(surfnr, new mesh.Point(pmap.Get(1)), pgi.Elem(1), n);
        
				for (int j = 1; j <= rule.oldels.Size(); j++)
				{
				  bad1 += mesh.SurfaceElement(elmap.Get(j)).CalcJacobianBadness(mesh.Points(), n);
				}
        
				// check new element:
				for (int j = 1; j <= rule.newels.Size(); j++)
				{
					Element2d rnel = rule.newels.Get(j);
					Element2d nel = new Element2d(rnel.GetNP());
					for (int k = 1; k <= rnel.GetNP(); k++)
					{
					  nel.PNum(k) = pmap.Get(rnel.PNum(k));
					}
        
					bad2 += nel.CalcJacobianBadness(mesh.Points(), n);
				}
        
				if (bad2 > 1e3)
				{
					continue;
				}
        
				if (newdef == olddef && bad2 > bad1)
				{
					continue;
				}
        
        
				// generate new element:
				for (int j = 1; j <= rule.newels.Size(); j++)
				{
					Element2d rnel = rule.newels.Get(j);
					Element2d nel = new Element2d(rnel.GetNP());
					nel.SetIndex(faceindex);
					for (int k = 1; k <= rnel.GetNP(); k++)
					{
					nel.PNum(k) = pmap.Get(rnel.PNum(k));
					nel.GeomInfoPi(k) = pgi.Get(rnel.PNum(k));
					}
        
					mesh.AddSurfaceElement(nel);
				}
        
				for (int j = 0; j < rule.oldels.Size(); j++)
				{
				  mesh.DeleteSurfaceElement(elmap[j]);
				}
        
				for (int j = 1; j <= pmap.Size(); j++)
				{
				  nelonnode[pmap.Get(j)] += rule.incelsonnode.Get(j);
				}
        
				used[ri]++;
				}
			}
			}
        
			mesh.Compress();
        
			for (int ri = 0; ri < rules.Size(); ri++)
			{
			PrintMessage(5, "rule ", ri + 1, " ", mapped[ri], "/", used[ri], " mapped/used");
			}
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void ProjectBoundaryPoints(Array<int> surfaceindex, Array<Point < 3> > from, Array<Point < 3> > dest)
		  {
			for (int i = 0; i < surfaceindex.Size(); i++)
			{
			if (surfaceindex[i] >= 0)
			{
				*dest[i] = *from[i];
				ProjectPoint(surfaceindex[i],*dest[i]);
			}
			}
        
        
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void ImproveVolumeMesh(Mesh mesh)
		  {
        
			if (!faceindex)
			{
			PrintMessage(3, "Smoothing");
        
			for (faceindex = 1; faceindex <= mesh.GetNFD(); faceindex++)
			{
				ImproveVolumeMesh(mesh);
				if (multithread.terminate)
				{
				  throw new Exception("Meshing stopped");
				}
			}
			faceindex = 0;
			return;
			}
        
        
        
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static int timer = NgProfiler::CreateTimer("MeshSmoothing 2D");
			NgProfiler.RegionTimer reg = new NgProfiler.RegionTimer(ImproveVolumeMesh_timer);
        
        
        
			CheckMeshApproximation(mesh);
        
			int i;
			int j;
			int k;
			SurfaceElementIndex sei = new SurfaceElementIndex();
        
			Array<SurfaceElementIndex> seia = new Array<SurfaceElementIndex>();
			mesh.GetSurfaceElementsOfFace(faceindex, seia);
        
			/*
			bool mixed = 0;
			for (i = 0; i < seia.Size(); i++)
			  if (mesh[seia[i]].GetNP() != 3)
			{
			  mixed = 1;
			  break;
			}
			*/
        
			int loci;
			double fact;
			bool moveisok;
        
			PointGeomInfo ngi = new PointGeomInfo();
			Point < 3> origp;
        
			Vector x = new Vector(3);
        
			Array<MeshPoint, PointIndex.BASE> savepoints = new Array<MeshPoint, PointIndex.BASE>(mesh.GetNP());
        
			Array<int, PointIndex.BASE> nelementsonpoint = new Array<int, PointIndex.BASE>(mesh.GetNP());
			nelementsonpoint = 0;
        
			for (i = 0; i < seia.Size(); i++)
			{
			Element2d el = mesh[seia[i]];
			for (j = 0; j < el.GetNP(); j++)
			{
			  nelementsonpoint[el[j]]++;
			}
			}
        
        
			TABLE<SurfaceElementIndex,PointIndex.BASE> elementsonpoint = new TABLE<SurfaceElementIndex,PointIndex.BASE>(nelementsonpoint);
			for (i = 0; i < seia.Size(); i++)
			{
			Element2d el = mesh[seia[i]];
			for (j = 0; j < el.GetNP(); j++)
			{
			  elementsonpoint.Add(el[j], seia[i]);
			}
			}
        
        
			JacobianPointFunction pf = new JacobianPointFunction(mesh.Points(), mesh.VolumeElements());
        
        
        
		//     Opti2SurfaceMinFunction surfminf(mesh);
		//     Opti2EdgeMinFunction edgeminf(mesh);
		//     Opti2SurfaceMinFunctionJacobian surfminfj(mesh);
        
			OptiParameters par = new OptiParameters();
			par.maxit_linsearch = 8;
			par.maxit_bfgs = 5;
        
			int np = mesh.GetNP();
			int ne = mesh.GetNE();
        
			BitArray badnodes = new BitArray(np);
			badnodes.Clear();
        
			for (i = 1; i <= ne; i++)
			{
			Element el = mesh.VolumeElement(i);
			double bad = el.CalcJacobianBadness(mesh.Points());
			if (bad > 1)
			{
			  for (j = 1; j <= el.GetNP(); j++)
			  {
				badnodes.Set(el.PNum(j));
			  }
			}
			}
        
        
			bool printeddot = false;
			char plotchar = '.';
			int modplot = 1;
			if (mesh.GetNP() > 1000)
			{
			plotchar = '+';
			modplot = 10;
			}
			if (mesh.GetNP() > 10000)
			{
			plotchar = 'o';
			modplot = 100;
			}
			int cnt = 0;
        
        
			Array<SurfaceElementIndex> locelements = new Array<SurfaceElementIndex>(0);
			Array<int> locrots = new Array<int>(0);
        
			for (PointIndex pi = mesh.Points().Begin(); pi < mesh.Points().End(); pi++)
			{
			if (mesh[pi].Type() != POINTTYPE.SURFACEPOINT)
			{
			  continue;
			}
        
			if (multithread.terminate)
			{
			  throw new Exception("Meshing stopped");
			}
        
			int surfi = -1;
        
			if (elementsonpoint[pi].Size() == 0)
			{
			  continue;
			}
        
			Element2d hel = mesh[elementsonpoint[pi][0]];
        
			if (hel.GetIndex() != faceindex)
			{
			  continue;
			}
        
			cnt++;
			if (cnt % modplot == 0 && writestatus)
			{
				printeddot = true;
				PrintDot(plotchar);
			}
        
        
			int hpi = 0;
			for (j = 1; j <= hel.GetNP(); j++)
			{
			  if (hel.PNum(j) == pi)
			  {
				  hpi = j;
				  break;
			  }
			}
			PointGeomInfo gi1 = hel.GeomInfoPi(hpi);
        
			locelements.SetSize(0);
			locrots.SetSize(0);
        
			for (j = 0; j < elementsonpoint[pi].Size(); j++)
			{
				sei = elementsonpoint[pi][j];
				Element2d bel = mesh[sei];
				surfi = mesh.GetFaceDescriptor(bel.GetIndex()).SurfNr();
        
				locelements.Append(sei);
        
				for (k = 1; k <= bel.GetNP(); k++)
				{
				  if (bel.PNum(k) == pi)
				  {
				  locrots.Append(k);
				  break;
				  }
				}
			}
        
        
			double lh = mesh.GetH(new mesh.Point(pi));
			par.typx = lh;
        
			pf.SetPointIndex(new netgen.PointIndex(pi));
        
			x = 0;
			bool pok = (pf.Func(x) < 1e10);
        
			if (pok)
			{
				BFGS(x, pf, par);
        
				origp = mesh[pi];
				loci = 1;
				fact = 1;
				moveisok = false;
        
        
				//optimizer loop (if whole distance is not possible, move only a bit!!!!)
				while (loci <= 5 && !moveisok)
				{
				loci++;
				mesh[pi](0) = origp(0) + x(0) * fact;
				mesh[pi](1) = origp(1) + x(1) * fact;
				mesh[pi](2) = origp(2) + x(2) * fact;
				fact = fact / 2.0;
        
        
				//cout << "origp " << origp << " newp " << mesh[pi];
        
				ngi.CopyFrom(gi1);
				moveisok = (ProjectPointGI(surfi, mesh[pi], ngi) != 0);
        
				//cout << " projected " << mesh[pi] << endl;
        
				// point lies on same chart in stlsurface
        
				if (moveisok)
				{
					for (j = 0; j < locelements.Size(); j++)
					{
					  mesh[locelements[j]].GeomInfoPi(locrots[j]) = ngi;
					}
        
					//cout << "moved " << origp << " to " << mesh[pi] << endl;
				}
				else
				{
					mesh[pi] = origp;
				}
        
				}
			}
			else
			{
				Console.Write("el not ok (point ");
				Console.Write(pi);
				Console.Write(": ");
				Console.Write(mesh[pi]);
				Console.Write(")");
				Console.Write("\n");
			}
			}
        
			if (printeddot)
			{
			  PrintDot('\n');
			}
        
			CheckMeshApproximation(mesh);
			mesh.SetNextTimeStamp();
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public MeshOptimize2d()
		  {
			SetFaceIndex(0);
			SetImproveEdges(0);
			SetMetricWeight(0);
			SetWriteStatus(1);
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void SelectSurfaceOfPoint(Point < 3> p, PointGeomInfo gi)
		  {
			;
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void ImproveMesh(Mesh mesh, MeshingParameters mp)
		  {
			if (!faceindex)
			{
			PrintMessage(3, "Smoothing");
        
			for (faceindex = 1; faceindex <= mesh.GetNFD(); faceindex++)
			{
				ImproveMesh(mesh, mp);
				if (multithread.terminate)
				{
				  throw new Exception("Meshing stopped");
				}
			}
			faceindex = 0;
			return;
			}
        
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static int timer = NgProfiler::CreateTimer("MeshSmoothing 2D");
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static int timer1 = NgProfiler::CreateTimer("MeshSmoothing 2D start");
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static int timer2 = NgProfiler::CreateTimer("MeshSmoothing 2D - BFGS");
        
			NgProfiler.RegionTimer reg = new NgProfiler.RegionTimer(ImproveMesh_timer);
			NgProfiler.StartTimer(ImproveMesh_timer1);
        
			CheckMeshApproximation(mesh);
        
			Opti2dLocalData ld = new Opti2dLocalData();
        
        
			Array<SurfaceElementIndex> seia = new Array<SurfaceElementIndex>();
			mesh.GetSurfaceElementsOfFace(faceindex, seia);
			bool mixed = false;
			for (int i = 0; i < seia.Size(); i++)
			{
			  if (mesh[seia[i]].GetNP() != 3)
			  {
			  mixed = true;
			  break;
			  }
			}
        
			Vector x = new Vector(2);
        
			Array<MeshPoint, PointIndex.BASE> savepoints = new Array<MeshPoint, PointIndex.BASE>(mesh.GetNP());
        
			ld.uselocalh = mp.uselocalh;
        
			Array<int, PointIndex.BASE> compress = new Array<int, PointIndex.BASE>(mesh.GetNP());
			Array<PointIndex> icompress = new Array<PointIndex>();
			for (int i = 0; i < seia.Size(); i++)
			{
			Element2d el = mesh[seia[i]];
			for (int j = 0; j < el.GetNP(); j++)
			{
			  compress[el[j]] = -1;
			}
			}
			for (int i = 0; i < seia.Size(); i++)
			{
			Element2d el = mesh[seia[i]];
			for (int j = 0; j < el.GetNP(); j++)
			{
			  if (compress[el[j]] == -1)
			  {
				  compress[el[j]] = icompress.Size();
				  icompress.Append(el[j]);
			  }
			}
			}
			Array<int> cnta = new Array<int>(icompress.Size());
			cnta = 0;
			for (int i = 0; i < seia.Size(); i++)
			{
			Element2d el = mesh[seia[i]];
			for (int j = 0; j < el.GetNP(); j++)
			{
			  cnta[compress[el[j]]]++;
			}
			}
			TABLE<SurfaceElementIndex> elementsonpoint = new TABLE<SurfaceElementIndex>(cnta);
			for (int i = 0; i < seia.Size(); i++)
			{
			Element2d el = mesh[seia[i]];
			for (int j = 0; j < el.GetNP(); j++)
			{
			  elementsonpoint.Add(compress[el[j]], seia[i]);
			}
			}
        
        
			/*
			Array<int, PointIndex::BASE> nelementsonpoint(mesh.GetNP());
			nelementsonpoint = 0;
			for (int i = 0; i < seia.Size(); i++)
			  {
			const Element2d & el = mesh[seia[i]];
			for (int j = 0; j < el.GetNP(); j++)
			  nelementsonpoint[el[j]]++;
			  }
		
			TABLE<SurfaceElementIndex,PointIndex::BASE> elementsonpoint(nelementsonpoint);
		
			for (int i = 0; i < seia.Size(); i++)
			  {
			const Element2d & el = mesh[seia[i]];
			for (int j = 0; j < el.GetNP(); j++)
			  elementsonpoint.Add (el[j], seia[i]);
			  }
			*/
        
        
        
			ld.loch = mp.maxh;
			ld.locmetricweight = metricweight;
			ld.meshthis = this;
        
        
        
			Opti2SurfaceMinFunction surfminf = new Opti2SurfaceMinFunction(mesh, ld);
			Opti2EdgeMinFunction edgeminf = new Opti2EdgeMinFunction(mesh, ld);
			Opti2SurfaceMinFunctionJacobian surfminfj = new Opti2SurfaceMinFunctionJacobian(mesh, ld);
        
			OptiParameters par = new OptiParameters();
			par.maxit_linsearch = 8;
			par.maxit_bfgs = 5;
        
			/*
			int i, j, k;
			Vector xedge(1);
			  if (improveedges)
			  for (i = 1; i <= mesh.GetNP(); i++)
			  if (mesh.PointType(i) == EDGEPOINT)
			  {
			  continue;
			  PrintDot ();
			  sp1 = mesh.Point(i);
		
			  locelements.SetSize(0);
			  locrots.SetSize (0);
			  lochs.SetSize (0);
			  surfi = surfi2 = surfi3 = 0;
		
			  for (j = 0; j < elementsonpoint[i].Size(); j++)
			  {
			  sei = elementsonpoint[i][j];
			  const Element2d * bel = &mesh[sei];
		
			  if (!surfi)
			  surfi = mesh.GetFaceDescriptor(bel->GetIndex()).SurfNr();
			  else if (surfi != mesh.GetFaceDescriptor(bel->GetIndex()).SurfNr())
			  {
			  if (surfi2 != 0 && surfi2 != 
			  mesh.GetFaceDescriptor(bel->GetIndex()).SurfNr())
			  surfi3 = mesh.GetFaceDescriptor(bel->GetIndex()).SurfNr();
			  else
			  surfi2 = mesh.GetFaceDescriptor(bel->GetIndex()).SurfNr();
			  }
		
			  locelements.Append (sei);
		
			  if (bel->PNum(1) == i)
			  locrots.Append (1);
			  else if (bel->PNum(2) == i)
			  locrots.Append (2);
			  else
			  locrots.Append (3);
		
			  if (uselocalh)
			  {
			  Point3d pmid = Center (mesh.Point(bel->PNum(1)),
			  mesh.Point(bel->PNum(2)),
			  mesh.Point(bel->PNum(3)));
			  lochs.Append (mesh.GetH(pmid));
			  }
			  }
		
			  if (surfi2 && !surfi3)
			  {
			  Vec3d n1, n2;
			  GetNormalVector (surfi, sp1, n1);
			  GetNormalVector (surfi2, sp1, n2);
			  t1 = Cross (n1, n2);
		
			  xedge = 0;
			  BFGS (xedge, edgeminf, par, 1e-6);
		
			  mesh.Point(i).X() += xedge.Get(1) * t1.X();
			  mesh.Point(i).Y() += xedge.Get(1) * t1.Y();
			  mesh.Point(i).Z() += xedge.Get(1) * t1.Z();
			  ProjectPoint2 (surfi, surfi2, mesh.Point(i));
			  }
			  }
			*/
        
        
			bool printeddot = false;
			char plotchar = '.';
			int modplot = 1;
			if (mesh.GetNP() > 1000)
			{
			plotchar = '+';
			modplot = 100;
			}
			if (mesh.GetNP() > 10000)
			{
			plotchar = 'o';
			modplot = 1000;
			}
			if (mesh.GetNP() > 100000)
			{
			plotchar = 'O';
			modplot = 10000;
			}
			int cnt = 0;
        
        
			NgProfiler.StopTimer(ImproveMesh_timer1);
        
			/*
			for (PointIndex pi = PointIndex::BASE; pi < mesh.GetNP()+PointIndex::BASE; pi++)
			  if (mesh[pi].Type() == SURFACEPOINT)
			*/
			for (int hi = 0; hi < icompress.Size(); hi++)
			{
			PointIndex pi = icompress[hi];
			if (mesh[pi].Type() == POINTTYPE.SURFACEPOINT)
			{
				if (multithread.terminate)
				{
				  throw new Exception("Meshing stopped");
				}
        
				cnt++;
				if (cnt % modplot == 0 && writestatus)
				{
				printeddot = true;
				PrintDot(plotchar);
				}
        
				// if (elementsonpoint[pi].Size() == 0) continue;
				if (elementsonpoint[hi].Size() == 0)
				{
					continue;
				}
        
				ld.sp1 = mesh[pi];
        
				// Element2d & hel = mesh[elementsonpoint[pi][0]];
				Element2d hel = mesh[elementsonpoint[hi][0]];
        
				int hpi = 0;
				for (int j = 1; j <= hel.GetNP(); j++)
				{
				  if (hel.PNum(j) == pi)
				  {
				  hpi = j;
				  break;
				  }
				}
        
				ld.gi1.CopyFrom(hel.GeomInfoPi(hpi));
				SelectSurfaceOfPoint(ld.sp1, ld.gi1);
        
				ld.locelements.SetSize(0);
				ld.locrots.SetSize(0);
				ld.lochs.SetSize(0);
					ld.loc_pnts2.SetSize(0);
					ld.loc_pnts3.SetSize(0);
        
				for (int j = 0; j < elementsonpoint[hi].Size(); j++)
				{
				SurfaceElementIndex sei = elementsonpoint[hi][j];
				Element2d bel = mesh[sei];
				ld.surfi = mesh.GetFaceDescriptor(bel.GetIndex()).SurfNr();
        
				ld.locelements.Append(sei);
        
				for (int k = 1; k <= bel.GetNP(); k++)
				{
				  if (bel.PNum(k) == pi)
				  {
					  ld.locrots.Append(k);
							  ld.loc_pnts2.Append(mesh[bel.PNumMod(k + 1)]);
							  ld.loc_pnts3.Append(mesh[bel.PNumMod(k + 2)]);
					  break;
				  }
				}
        
				if (ld.uselocalh != 0)
				{
					Point3d pmid = Center(mesh[bel[0]], mesh[bel[1]], mesh[bel[2]]);
					ld.lochs.Append(mesh.GetH(pmid));
				}
				}
        
			  GetNormalVector(ld.surfi, ld.sp1, ld.gi1, ld.normal);
			  ld.t1 = ld.normal.GetNormal();
			  ld.t2 = Cross(ld.normal, ld.t1.functorMethod);
        
			  // save points, and project to tangential plane
			  for (int j = 0; j < ld.locelements.Size(); j++)
			  {
				  Element2d el = mesh[ld.locelements[j]];
				  for (int k = 0; k < el.GetNP(); k++)
				  {
				savepoints[el[k]] = mesh[el[k]];
				  }
			  }
        
			  for (int j = 0; j < ld.locelements.Size(); j++)
			  {
				  Element2d el = mesh[ld.locelements[j]];
				  for (int k = 0; k < el.GetNP(); k++)
				  {
				  PointIndex hhpi = el[k];
				  double lam = ld.normal * (mesh[hhpi] - ld.sp1);
				  mesh[hhpi] -= lam * ld.normal;
				  }
			  }
        
			  x = 0;
			  par.typx = 0.3 * ld.lochs[0];
        
				  NgProfiler.StartTimer(ImproveMesh_timer2);
        
			  if (mixed)
			  {
				  BFGS(x, surfminfj, par, 1e-6);
			  }
			  else
			  {
				  BFGS(x, surfminf, par, 1e-6);
			  }
        
				  NgProfiler.StopTimer(ImproveMesh_timer2);
        
			  Point3d origp = mesh[pi];
			  int loci = 1;
			  double fact = 1;
			  int moveisok = 0;
        
			  // restore other points
			  for (int j = 0; j < ld.locelements.Size(); j++)
			  {
				  Element2d el = mesh[ld.locelements[j]];
				  for (int k = 0; k < el.GetNP(); k++)
				  {
				  PointIndex hhpi = el[k];
				  if (hhpi != pi)
				  {
					  mesh[hhpi] = savepoints[hhpi];
				  }
				  }
			  }
        
        
			  //optimizer loop (if whole distance is not possible, move only a bit!!!!)
			  while (loci <= 5 && moveisok == 0)
			  {
				  loci++;
					  /*
				  mesh[pi].X() = origp.X() + (x.Get(1) * t1.X() + x.Get(2) * t2.X())*fact;
				  mesh[pi].Y() = origp.Y() + (x.Get(1) * t1.Y() + x.Get(2) * t2.Y())*fact;
				  mesh[pi].Z() = origp.Z() + (x.Get(1) * t1.Z() + x.Get(2) * t2.Z())*fact;
					  */
					  Vec < 3> hv = x(0) * ld.t1.functorMethod + x(1) * ld.t2.functorMethod;
					  Point3d hnp = origp + new Vec3d(hv);
					  mesh[pi](0) = hnp.X();
					  mesh[pi](1) = hnp.Y();
					  mesh[pi](2) = hnp.Z();
        
				  fact = fact / 2.0;
        
				  // ProjectPoint (surfi, mesh[pi]);
				  // moveisok = CalcPointGeomInfo(surfi, ngi, mesh[pi]); 
        
				  PointGeomInfo ngi = new PointGeomInfo();
				  ngi.CopyFrom(ld.gi1);
				  moveisok = ProjectPointGI(ld.surfi, mesh[pi], ngi);
				  // point lies on same chart in stlsurface
        
				  if (moveisok != 0)
				  {
				  for (int j = 0; j < ld.locelements.Size(); j++)
				  {
					mesh[ld.locelements[j]].GeomInfoPi(ld.locrots[j]) = ngi;
				  }
				  }
				  else
				  {
				  mesh[pi] = Point < 3> (origp);
				  }
        
			  }
			}
			}
			if (printeddot)
			{
			  PrintDot('\n');
			}
        
			CheckMeshApproximation(mesh);
			mesh.SetNextTimeStamp();
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void GetNormalVector(int UnnamedParameter, Point < 3> p, ref Vec < 3> nv)
		  {
			nv = Vec < 3> (0, 0, 1);
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void GetNormalVector(int surfind, Point < 3> p, PointGeomInfo gi, Vec < 3> n)
		  {
			GetNormalVector(surfind, p.functorMethod, n.functorMethod);
		  }
	}
}