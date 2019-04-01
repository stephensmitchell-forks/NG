namespace netgen
{

	public class LocalH
	{
//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void FindInnerBoxes(AdFront3 adfront, testinnerDelegate testinner)
		  {
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static int timer = NgProfiler::CreateTimer("LocalH::FindInnerBoxes");
			NgProfiler.RegionTimer reg = new NgProfiler.RegionTimer(FindInnerBoxes_timer);
        
        
			int nf = adfront.GetNF();
        
			for (int i = 0; i < boxes.Size(); i++)
			{
			  boxes[i].flags.isinner = 0;
			}
        
			root.flags.isinner = 0;
        
			Point3d rpmid = new Point3d(root.xmid[0], root.xmid[1], root.xmid[2]);
			Vec3d rv = new Vec3d(root.h2, root.h2, root.h2);
			Point3d rx2 = rpmid + rv;
			// Point3d rx1 = rpmid - rv;
        
        
			root.flags.pinner = !adfront.SameSide(rpmid, rx2);
        
			if (testinner != null)
			{
			  (*testout) << "inner = " << root.flags.pinner << " =?= " << testinner(new Point3d(root.xmid[0], root.xmid[1], root.xmid[2])) << "\n";
			}
        
			Array<int> faceinds = new Array<int>(nf);
			Array<Box3d> faceboxes = new Array<Box3d>(nf);
        
			for (int i = 1; i <= nf; i++)
			{
			faceinds.Elem(i) = i;
			adfront.GetFaceBoundingBox(i, faceboxes.Elem(i));
			}
        
			for (int i = 0; i < 8; i++)
			{
			  FindInnerBoxesRec2(root.childs[i], adfront, faceboxes, faceinds, nf);
			}
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void FindInnerBoxes(AdFront2 adfront, testinnerDelegate testinner)
		  {
		//C++ TO C# CONVERTER NOTE: This static local variable declaration (not allowed in C#) has been moved just prior to the method:
		//	static int timer = NgProfiler::CreateTimer("LocalH::FindInnerBoxes 2d");
			NgProfiler.RegionTimer reg = new NgProfiler.RegionTimer(FindInnerBoxes_timer);
        
			for (int i = 0; i < boxes.Size(); i++)
			{
			  boxes[i].flags.isinner = 0;
			}
        
			root.flags.isinner = 0;
        
			Point < 2> rpmid(root.xmid[0], root.xmid[1]); // , root->xmid[2]);
			Vec < 2> rv(root.h2, root.h2);
			Point < 2> rx2 = rpmid + rv;
			// Point<2> rx1 = rpmid - rv;
        
        
			root.flags.pinner = !adfront.SameSide(rpmid, rx2);
        
			if (testinner != null)
			{
			  (*testout) << "inner = " << root.flags.pinner << " =?= " << testinner(rpmid) << "\n";
			}
        
        
			int nf = adfront.GetNFL();
			Array<int> faceinds = new Array<int>(nf);
			Array<Box < 3>> faceboxes = new Array<Box < 3>>(nf);
        
			for (int i = 0; i < nf; i++)
			{
			faceinds[i] = i;
			// adfront->GetFaceBoundingBox(i, faceboxes.Elem(i));
        
			FrontLine line = adfront.GetLine(i);
			faceboxes[i].Set(adfront.GetPoint(line.L().I1()));
			faceboxes[i].Add(adfront.GetPoint(line.L().I2()));
			}
        
			for (int i = 0; i < 8; i++)
			{
			  FindInnerBoxesRec2(root.childs[i], adfront, faceboxes, faceinds, nf);
			}
		  }
	}
}