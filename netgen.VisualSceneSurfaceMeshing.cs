namespace netgen
{

	public class VisualSceneSurfaceMeshing
	{
//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public VisualSceneSurfaceMeshing()
		  {
			  this.VisualScene = new <type missing>();
			;
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void Dispose()
		  {
			;
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void DrawScene()
		  {
			int i;
			int j;
			int k;
        
			if (loclines.Size() != changeval)
			{
			center = Point < 3>(0,0,-5);
			rad = 0.1;
        
			CalcTransformationMatrices();
			changeval = loclines.Size();
			}
        
		  glClearColor(backcolor, backcolor, backcolor, 1.0);
		  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		  SetLight();
        
		  //  glEnable (GL_COLOR_MATERIAL);
        
		  //  glDisable (GL_SHADING);
		  //  glColor3f (0.0f, 1.0f, 1.0f);
		  //  glLineWidth (1.0f);
		  //  glShadeModel (GL_SMOOTH);
        
		  //  glCallList (linelists.Get(1));
        
		  //  SetLight();
        
		  glPushMatrix();
		  glMultMatrixf(transformationmat);
        
		  glShadeModel(GL_SMOOTH);
		  // glDisable (GL_COLOR_MATERIAL);
		  glEnable(GL_COLOR_MATERIAL);
		  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        
		  glEnable(GL_BLEND);
		  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        
		  //  glEnable (GL_LIGHTING);
        
		  double shine = vispar.shininess;
		  double transp = vispar.transp;
        
		  glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, shine);
		  glLogicOp(GL_COPY);
        
        
        
		  /*
		
		  float mat_col[] = { 0.2, 0.2, 0.8, 1 };
		  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_col);
		
		  glPolygonOffset (1, 1);
		  glEnable (GL_POLYGON_OFFSET_FILL);
		
		    float mat_colbl[] = { 0.8, 0.2, 0.2, 1 };
		    float mat_cololdl[] = { 0.2, 0.8, 0.2, 1 };
		    float mat_colnewl[] = { 0.8, 0.8, 0.2, 1 };
		
		
		    glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
		    glPolygonOffset (1, -1);
		    glLineWidth (3);
		
		    for (i = 1; i <= loclines.Size(); i++)
		      {
			if (i == 1)
			  {
				glEnable (GL_POLYGON_OFFSET_FILL);
				glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_colbl);
			  }
			else if (i <= oldnl) 
			  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_cololdl);
			else 
			  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_colnewl);
		
			int pi1 = loclines.Get(i).I1();
			int pi2 = loclines.Get(i).I2();
		
			if (pi1 >= 1 && pi2 >= 1)
			  {
				Point3d p1 = locpoints.Get(pi1);
				Point3d p2 = locpoints.Get(pi2);
		
				glBegin (GL_LINES);
				glVertex3f (p1.X(), p1.Y(), p1.Z());
				glVertex3f (p2.X(), p2.Y(), p2.Z());
				glEnd();
			  }
		
			glDisable (GL_POLYGON_OFFSET_FILL);
		      }
		  
		
		    glLineWidth (1);
		
		
		    glPointSize (5);
		    float mat_colp[] = { 1, 0, 0, 1 };
		    glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_colp);
		    glBegin (GL_POINTS);
		    for (i = 1; i <= locpoints.Size(); i++)
		      {
			Point3d p = locpoints.Get(i);
			glVertex3f (p.X(), p.Y(), p.Z());
		      }
		    glEnd();
		
		
		    glPopMatrix();
		  */
        
			float[] mat_colp = {1F, 0F, 0F, 1F};
        
			float[] mat_col2d1 = {1F, 0.5F, 0.5F, 1F};
			float[] mat_col2d = {1F, 1F, 1F, 1F};
			glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_col2d);
        
			double scalex = 0.1;
			double scaley = 0.1;
        
			glBegin(GL_LINES);
			for (i = 1; i <= loclines.Size(); i++)
			{
			glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_col2d);
			if (i == 1)
			{
			  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_col2d1);
			}
        
			int pi1 = loclines.Get(i).I1();
			int pi2 = loclines.Get(i).I2();
        
			if (pi1 >= 1 && pi2 >= 1)
			{
				Point2d p1 = plainpoints.Get(pi1);
				Point2d p2 = plainpoints.Get(pi2);
        
				glBegin(GL_LINES);
				glVertex3f(scalex * p1.X(), scaley * p1.Y(), -5);
				glVertex3f(scalex * p2.X(), scaley * p2.Y(), -5);
				glEnd();
			}
			}
			glEnd();
        
        
			glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_colp);
			glBegin(GL_POINTS);
			for (i = 1; i <= plainpoints.Size(); i++)
			{
			Point2d p = plainpoints.Get(i);
			glVertex3f(scalex * p.X(), scaley * p.Y(), -5);
			}
			glEnd();
        
        
        
        
        
        
		  glDisable(GL_POLYGON_OFFSET_FILL);
        
		  glPopMatrix();
		  DrawCoordinateCross();
		  DrawNetgenLogo();
		  glFinish();
        
		  /*
		    glDisable (GL_POLYGON_OFFSET_FILL);
		
		    //  cout << "draw surfacemeshing" << endl;
		    //
		    //  if (changeval != stlgeometry->GetNT())
		    //      BuildScene();
		    //      changeval = stlgeometry->GetNT();
		    
		
		    glClearColor(backcolor, backcolor, backcolor, 1.0);
		    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		
		    SetLight();
		
		    glPushMatrix();
		    glLoadMatrixf (transmat);
		    glMultMatrixf (rotmat);
		
		    glShadeModel (GL_SMOOTH);
		    glDisable (GL_COLOR_MATERIAL);
		    glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
		
		    glEnable (GL_BLEND);
		    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		
		    float mat_spec_col[] = { 1, 1, 1, 1 };
		    glMaterialfv (GL_FRONT_AND_BACK, GL_SPECULAR, mat_spec_col);
		
		    double shine = vispar.shininess;
		    double transp = vispar.transp;
		
		    glMaterialf (GL_FRONT_AND_BACK, GL_SHININESS, shine);
		    glLogicOp (GL_COPY);
		
		
		    float mat_col[] = { 0.2, 0.2, 0.8, transp };
		    float mat_colrt[] = { 0.2, 0.8, 0.8, transp };
		    glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_col);
		
		    glPolygonOffset (1, 1);
		    glEnable (GL_POLYGON_OFFSET_FILL);
		
		    glColor3f (1.0f, 1.0f, 1.0f);
		
		    glEnable (GL_NORMALIZE);
		    
		    //  glBegin (GL_TRIANGLES);
		    //      for (j = 1; j <= stlgeometry -> GetNT(); j++)
		    //      {
		    //      glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_col);
		    //      if (j == geomtrig)
		    //      glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_colrt);
		
		
		    //      const STLReadTriangle & tria = stlgeometry -> GetReadTriangle(j);
		    //      glNormal3f (tria.normal.X(),
		    //      tria.normal.Y(),
		    //      tria.normal.Z());
		
		    //      for (k = 0; k < 3; k++)
		    //      {
		    //      glVertex3f (tria.pts[k].X(),
		    //      tria.pts[k].Y(),
		    //      tria.pts[k].Z());
		    //      }
		    //      }    
		    //      glEnd ();
		    
		
		
		    glDisable (GL_POLYGON_OFFSET_FILL);
		
		    float mat_colbl[] = { 0.8, 0.2, 0.2, 1 };
		    float mat_cololdl[] = { 0.2, 0.8, 0.2, 1 };
		    float mat_colnewl[] = { 0.8, 0.8, 0.2, 1 };
		
		
		    glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
		    glPolygonOffset (1, -1);
		    glLineWidth (3);
		
		    for (i = 1; i <= loclines.Size(); i++)
		      {
			if (i == 1)
			  {
				glEnable (GL_POLYGON_OFFSET_FILL);
				glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_colbl);
			  }
			else if (i <= oldnl) 
			  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_cololdl);
			else 
			  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_colnewl);
		
			int pi1 = loclines.Get(i).I1();
			int pi2 = loclines.Get(i).I2();
		
			if (pi1 >= 1 && pi2 >= 1)
			  {
				Point3d p1 = locpoints.Get(pi1);
				Point3d p2 = locpoints.Get(pi2);
		
				glBegin (GL_LINES);
				glVertex3f (p1.X(), p1.Y(), p1.Z());
				glVertex3f (p2.X(), p2.Y(), p2.Z());
				glEnd();
			  }
		
			glDisable (GL_POLYGON_OFFSET_FILL);
		      }
		
		
		    glLineWidth (1);
		
		
		    glPointSize (5);
		    float mat_colp[] = { 1, 0, 0, 1 };
		    glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_colp);
		    glBegin (GL_POINTS);
		    for (i = 1; i <= locpoints.Size(); i++)
		      {
			Point3d p = locpoints.Get(i);
			glVertex3f (p.X(), p.Y(), p.Z());
		      }
		    glEnd();
		
		
		    glPopMatrix();
		
		
		    float mat_col2d1[] = { 1, 0.5, 0.5, 1 };
		    float mat_col2d[] = { 1, 1, 1, 1 };
		    glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_col2d);
		  
		    double scalex = 0.1, scaley = 0.1;
		
		    glBegin (GL_LINES);
		    for (i = 1; i <= loclines.Size(); i++)
		      {
			glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_col2d);
			if (i == 1)
			  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_col2d1);
		
			int pi1 = loclines.Get(i).I1();
			int pi2 = loclines.Get(i).I2();
		
			if (pi1 >= 1 && pi2 >= 1)
			  {
				Point2d p1 = plainpoints.Get(pi1);
				Point2d p2 = plainpoints.Get(pi2);
		
				glBegin (GL_LINES);
				glVertex3f (scalex * p1.X(), scaley * p1.Y(), -5);
				glVertex3f (scalex * p2.X(), scaley * p2.Y(), -5);
				glEnd();
			  }
		      }
		    glEnd ();
		
		
		    glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_colp);
		    glBegin (GL_POINTS);
		    for (i = 1; i <= plainpoints.Size(); i++)
		      {
			Point2d p = plainpoints.Get(i);
			glVertex3f (scalex * p.X(), scaley * p.Y(), -5);
		      }
		    glEnd();
		
		    glFinish();  
		*/
		  }

//C++ TO C# CONVERTER WARNING: The original C++ declaration of the following method implementation was not found:
		  public void BuildScene(int zoomall)
		  {
			int i;
			int j;
			int k;
			/*
			  center = stlgeometry -> GetBoundingBox().Center();
			  rad = stlgeometry -> GetBoundingBox().Diam() / 2;
		
			  CalcTransformationMatrices();
			*/
		  }
	}
}