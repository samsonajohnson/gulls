#include<vector>
#include<iostream>
#include<fstream>
#include<cmath>
#include<fitsio.h>
#include<ctime>
#include<cstdlib>
#include<string>
#include<sstream>
// #include<algorithm>   //added initially for min/max elements in course vectors
// #include<string>   //added for double to string conversion in file naming

#include "psf.h"
#include "image.h"
#include "constants.h"
#include "split.h"
#include "integerPowers.h"

// Variables for use with course data
// Easier to define globally and read for each iteration of images, even if it's not
// the most efficient
vector<double> course_year, course_x, course_y, course_theta;
int bounds_xmin, bounds_xmax, bounds_ymin, bounds_ymax;
bool bounds_set = false;
double origin_year, origin_x, origin_y, origin_theta;


double fs(double Fin, double Fout)
{
  double Fsource = Fin-Fout;
  return Fsource/Fin;
}

void load_course(string file_name)
{
  // clear course data
  course_year.clear();
  course_x.clear();
  course_y.clear();
  course_theta.clear();

  //load dithering course file
  ifstream coursein(file_name.c_str());
  if(!coursein)
  {
    cout << "Could not open dithering course file: " << file_name << endl;
    exit(1);
  }

  while(!coursein.eof())
  {
    string line;
    vector<string> data;

    getline(coursein,line);
    split(line,data);

    if(int(data.size())>=4)
    {
      course_year.push_back(atof(data[0].c_str())/365.25);
      course_x.push_back(atof(data[1].c_str()));
      course_y.push_back(atof(data[2].c_str()));
      course_theta.push_back(atof(data[3].c_str()));
    }
  }

  coursein.close();

  //normalizes course_x and course_y so that first waypoint is treated as 0,0
  origin_year = course_year[0];
  origin_x = course_x[0];
  origin_y = course_y[0];
  origin_theta = course_theta[0];

  for(int i=0;i<course_x.size();i++)
  {
    course_year[i] -= origin_year;
    course_x[i] -= origin_x;
    course_y[i] -= origin_y;
    course_theta[i] -= origin_theta;
  }

  int course_xmin = int(floor(*min_element(course_x.begin(),course_x.end())));
  int course_xmax = int(ceil(*max_element(course_x.begin(),course_x.end())));
  int course_ymin = int(floor(*min_element(course_y.begin(),course_y.end())));
  int course_ymax = int(ceil(*max_element(course_y.begin(),course_y.end())));

  if(course_xmin<bounds_xmin)
  {
    bounds_xmin = course_xmin;
  }

  if(course_xmax>bounds_xmax)
  {
    bounds_xmax = course_xmax;
  }

  if(course_ymin<bounds_ymin)
  {
    bounds_ymin = course_ymin;
  }

  if(course_ymax>bounds_ymax)
  {
    bounds_ymax = course_ymax;
  }
}

int main(int argc, char* argv[])
{
//////////////////////////////////////////////////////////////////////////////
//                               USER OPTIONS
//

  const int maximsize = 1024;            //maximum image size without 
                                         //receiving a warning requiring
                                         //user interaction
  const bool output_image = true;        //output simulated images
  const double large_psf_mag = 15.0;     //magnitude where large psf is used
                                         //missed flux still gets distributed, 
                                         //but uniformly
  const int subpix_method = 0;           //Method by which sub-pixel systematic
                                         //is computed (see below for options)
  const bool zero_systematic = true;     //Do not include a systematic 
                                         //component in the raw error
                                         //can then add either const or 
                                         //calculated systematic

  //Subpixel computation:
  // 0 - the fractional rms fluctuation due the total flux within a slightly 
  //     expanded aperture. (USE THIS)
  // 1 - the fractional rms fluctuation due to the total flux of all stars 
  //     that lie within the same expanded aperture (DO NOT USE THIS)
  //
  //  i.e. method 1 gives a larger error if there are bright stars at the edge 
  //       of the aperture. Method 1 requires significantly more computation
  //       as for each photometered star it must look through all other stars
  //       to see if they are within the aperture. Method 0 just counts up the
  //       total flux in the expanded aperture
                                 
//  
//                               USER OPTIONS
/////////////////////////////////////////////////////////////////////////////

  if(argc<6 || (argc>=6 && argc%2==1))
  {
    cerr << "\nUsage:\n./transitfield <detectorlist> <nfilters> <field> <area> {<field> <area> {...} } <output filename root>\n" << endl;
    exit(1);
  }

  string detectorList = string(argv[1]);      //detector parameter file
  string outroot = string(argv[argc-1]);  //root of the output filename
  int nfilters=atoi(argv[2]);      //number of filters in the galactic model files
   
  //Output consists of:
  ////                   <outroot>.txt -the list of stars and photometry
  ////                   <outroot>_true.fits -true starfield image
  ////                   <outroot>_image.fits -single simulated starfield image
  ////                   <outroot>_stack.fits -stacked sim starfield image

  vector<string> detectorNames;
  vector<double> texp;
  vector<int> nstack;
  vector<int> xpix,ypix;
  vector<int> filter;
  vector<string> dither;

  double fieldrot=30*pi/180.0;
  double sfr=sin(fieldrot);
  double cfr=cos(fieldrot);

  //load the detectors

  string line;
  vector<string> data;

  ifstream detin(detectorList.c_str());
  if(!detin)
  {
    cerr << "Could not open the list of detectors: " << detectorList << endl;
    exit(1);
  }

  double largestx, largesty;
  int largestfilter = 0;

  while(!detin.eof())
  {
    getline(detin,line);
    split(line,data);

    if(int(data.size())>=6)
    {
      detectorNames.push_back(data[0]);
      texp.push_back(atof(data[1].c_str()));
      nstack.push_back(atoi(data[2].c_str()));
      xpix.push_back(atoi(data[3].c_str()));
      ypix.push_back(atoi(data[4].c_str()));
      filter.push_back(atoi(data[5].c_str()));
      dither.push_back(data[6]);

      if(bounds_set==false)
      {
        //load dithering course file
        string file_name = dither[dither.size()-1];
        ifstream coursein(file_name.c_str());
        if(!coursein)
        {
          cout << "Could not open dithering course file: " << file_name << endl;
          exit(1);
        }

        string line;
        vector<string> data;

        getline(coursein,line);
        split(line,data);

        if(int(data.size())>=4)
        {
          bounds_xmin = atof(data[1].c_str());
          bounds_xmax = atof(data[1].c_str());
          bounds_ymin = atof(data[2].c_str());
          bounds_ymax = atof(data[2].c_str());
        }

        coursein.close();
        bounds_set = true;

      }
      load_course(dither[dither.size()-1]);   // For setting min/max bounds of dither area

      if(filter[filter.size()-1]>largestfilter) largestfilter = filter[filter.size()-1];
    }
  }

  string course_file, output_file;
  stringstream i_buff;
  for(int i=0; i<dither.size(); i++)
  {
    course_file = dither[i];
    ifstream src(course_file.c_str(), ios_base::binary);
    i_buff.clear();
    output_file = "";
    i_buff << "detector" << i << ".course";
    i_buff >> output_file;
    ofstream dst(output_file.c_str(), ios_base::binary);

    dst << src.rdbuf();
  }

  cout << detectorNames.size() << " observatories" << endl;

  vector<image> images(detectorNames.size());

  long seed;
  

  for(int im=0;im<int(detectorNames.size());im++)
  {
    if(images[im].load_detector(detectorNames[im])<0) exit(1);
    
    //free up the memory used by the psf loading process
    images[im].psf.flush_splines();

    //set the magnitude above which the small psf gets used
    images[im].set_largepsfmag(large_psf_mag);

    //set up the random seed
    images[im].pass_seed(&seed);
    if(im==0) images[im].reseed();

    images[im].set_image_properties(xpix[im],ypix[im]);

  }
  
  

  vector<string> fields; //filename of all the fields
  vector<double> areas;  //areas of all the fields
  vector<int> nstars;    //number of stars in all the fields
  vector<int> sid;    //pointer between star and its data
  vector<int> sfn;    //pointer between star and its data

  vector<int> paddock; //TOO MANY FIELDS!!!! The field number of each 
                       //photometered star
  vector<int> ref;     //reference to the star in the input list 

  int nfields=0;       //number of starfields that will get read in


  for(int i=3;i<argc-1;i+=2)
  {
    fields.push_back(string(argv[i]));
    areas.push_back(atof(argv[i+1]));
    nstars.push_back(0);
    nfields++;
  }

  //read in each starfield

  vector<vector<string> > stardata = vector<vector<string> >(nfields);

  //ignore any lines with these characters
  string ignore = 
    string("ABCDFGHIJKLMNOPQRSTUVWXYZabcdfghijklmnopqrstuvwxyz*:;=(),[]{}#");

  for(int i=0;i<nfields;i++)
  {
    ifstream in(fields[i].c_str());
    if(!in)
    {
      cerr << "Could not open field (" << fields[i] << ")." << endl;
      exit(1);
    }

    vector<double> data;
    string line;

    while(!in.eof())
    {
      getline(in,line);
      line = line.substr(0,line.find_first_of("#")); //remove comments

      if(line.find_first_of(ignore)==line.npos)
      {
        split(line,data);

        //keep track of the star's properties 
        //so we can do photometry with it later
        if(data.size()>=nfilters)
        {
          stardata[i].push_back(line);
          nstars[i]++;
        }
      }
    }

    cout << "Will use " << nstars[i] << " stars for photometry in field " << fields[i] << " over an area " << areas[i] << endl;

    in.close();
  }
  
  //start creating the image

  freestarlist stars;

  for(int i=0;i<int(detectorNames.size());i++) images[i].addbg();

  ofstream out((outroot + string("txt")).c_str());
  if(!out)
  {
    cerr << "Could not open output file (" << outroot + string(".txt") << ")" << endl;
    exit(1);
  }

  double scaledif;
  vector<double> psfback; //for each image
  int xmin,xmax,ymin,ymax;
  double mag;
  int setdumb, ndraw, j;   //vars for poisson draw process

  //Now add stars to the first image
  for(int i=0;i<nfields;i++)
  {
    areas[i]*=3600*3600; //convert to sq arcsec

    if(areas[i]<=0)
    {
      cerr << "Nonsense solid angle inputed (" << areas[i] << "). No stars were added to the image" << endl;
      exit(1);
    }

    //calculate the dimensions of the star-field required

    double Afield;

    cout << "field " << i << endl;

    load_course(dither[0]);
    xmin = bounds_xmin;
    xmax = bounds_xmax;
    ymin = bounds_ymin;
    ymax = bounds_ymax;

    images[0].field_dimensions(areas[i],xmin,xmax,ymin,ymax,Afield);
    
    cout << "stardata.size = " << stardata[i].size() << endl;

    //now we can generate the stars
    setdumb=0;
    ndraw = poisson(nstars[i]*Afield/areas[i],&seed);
    for(int reps=0; reps<ndraw; reps++)
	  {
      j = randint(0,nstars[i]-1,&seed);
      split(stardata[i][j],data);
      mag = atof(data[filter[0]].c_str());

      images[0].freeaddstar(xmin,xmax,ymin,ymax,mag,&stars);
      sfn.push_back(i);
      sid.push_back(j);

      for(int im=1;im<int(detectorNames.size());im++)
      {
        if(!setdumb)
        {
          cout << images[im].psf.Nsub << " " << images[0].psf.Nsub << endl;
        }
        // load_course(dither[im]);
        // xmin = bounds_xmin;
        // xmax = bounds_xmax;
        // ymin = bounds_ymin;
        // ymax = bounds_ymax;
        // images[im].field_dimensions(areas[i],xmin,xmax,ymin,ymax,Afield);

        scaledif=(images[im].psf.pixscale/images[0].psf.pixscale);
        mag = atof(data[filter[im]].c_str());

        //why is this needed? it's just one star
        //add star at a random position
        images[im].freeaddstar(stars.x[stars.nstars-1]/scaledif,stars.y[stars.nstars-1]/scaledif,mag);
      }
      setdumb=1;

      out << stars.x[stars.nstars-1] << " " << stars.y[stars.nstars-1] << " " << stardata[i][j] << "\n";
	  }

    //add the background due to missed psf tails
    for(int im=0;im<int(detectorNames.size());im++)
  	{
	    scaledif=(images[im].psf.pixscale/images[0].psf.pixscale);
	    //psfback[im] = images[im].addpsfbg(xmin/scaledif,xmax/scaledif,ymin/scaledif,ymax/scaledif);
	    images[im].addbg();
	  }

  }

  cout << "Added " << stars.nstars << " stars to the image" << endl;

  for(int im=0;im<int(detectorNames.size());im++)
  {
    char imno[10]; sprintf(imno,"%d",im);
    images[im].expose(texp[im],nstack[im]);
    images[im].write_fits(outroot+string("_detector")+string(imno)+string("_origin")+string("_stack.fits"),true);
  }



  double x0,y0;
  double mul, mub;
  double scale;
  double star_x, star_y;

  // stringstream ss;
  string dither_index;
  int init_len;
  
  //not sure why these have to be initialized to a particular value
  double trans_x=30*9;
  double trans_y=30*9;
  double x=30*9;
  double y=30*9;

  for(int im=0;im<int(detectorNames.size());im++)
  {
    cout << "  Detector " << im << endl;

    load_course(dither[im]);

    for(int w=0;w<course_x.size();w++)
    {
      cout << "    Waypoint " << w << endl;

      scale = images[im].psf.Nsub/(images[im].psf.pixscale*1000.0);
      images[im].reset_image();
      images[im].reset_detector();
      scaledif=(images[im].psf.pixscale/images[0].psf.pixscale);
      //y+=0.1372456;
      //images[im].freeaddstar(x,y,mag);

      double cos_theta = cos(course_theta[w]*pi/180.0);
      double sin_theta = sin(course_theta[w]*pi/180.0);

      for(int s=0;s<stars.nstars;s++)
      {
        split(stardata[sfn[s]][sid[s]],data);
        mag = atof(data[filter[0]].c_str());

        // scale converts from milliarcseconds/year to arcseconds/year
        // then we divide by 3600 to get degrees/year
        mul = atof(data[nfilters].c_str()) * scale;
        mub = atof(data[nfilters+1].c_str()) * scale;


        // Updating star locations for a dither position has 3 steps:

        // 1. Calculate location after proper motion
        star_x = stars.x[s] + course_year[w]*mul;
        star_y = stars.y[s] + course_year[w]*mub;

        // 2. Adjust for dither x,y offset
        trans_x = (star_x-course_x[w]*images[im].psf.Nsub)/scaledif;
        trans_y = (star_y-course_y[w]*images[im].psf.Nsub)/scaledif;

        // 3. Rotate around image center
        double x_cen = images[im].Xpix/2*images[im].psf.Nsub;
        double y_cen = images[im].Ypix/2*images[im].psf.Nsub;
        x = trans_x*cos_theta+trans_y*sin_theta-x_cen*cos_theta-y_cen*sin_theta+x_cen;
        y = trans_x*-1*sin_theta+trans_y*cos_theta+x_cen*sin_theta-y_cen*cos_theta+y_cen;


        images[im].freeaddstar(x, y, mag);
      }

      images[im].addbg();
      char imno[10]; sprintf(imno,"%d",im);
      
      char dithno[15]; sprintf(dithno,"%d",w);
      dither_index = string(dithno);
      // ss.str(string());
      // ss << w;
      // dither_index = "";
      // ss >> dither_index;
      init_len = dither_index.length();
      for(int sc=0; sc<6-init_len; sc++)
      {
        dither_index = "0" + dither_index;
      }

      images[im].expose(texp[im],nstack[im]);
      images[im].write_fits(outroot+string("_detector")+string(imno)+string("_")+dither_index+string("_stack.fits"),true);
      //images[im].freesubstar(x,y,mag);
    
    }
  }

  // cout << images[0].Xpix/2 << "   " << images[0].Ypix/2 << endl;

  cout << "Images written" << endl;

  // cout << "Course file had " << course_x.size() << " x entries and " << course_y.size() << " y entries" << endl;
  // for(int i=0;i<course_x.size();i++)
  // {
  //   cout << course_x[i] << " " << course_y[i] << endl;
  // }

}
