// Compile with:
//     
//
// To specify the number of bodies in the world, the program optionally accepts
// an integer as its first command line argument.

#include <time.h>
#include <sys/times.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <X11/Xlib.h>
#include <unistd.h>
#include <stdbool.h>
#include <omp.h>

#define WIDTH 1024
#define HEIGHT 768
#define THETA 0.0000005

// default number of bodies
#define DEF_NUM_BODIES 200
// gravitational constant
#define GRAV 10.0
// initial velocities are scaled by this value
#define V_SCALAR 20.0
// initial masses are scaled by this value
#define M_SCALAR 5.0
// radius scalar
#define R_SCALAR 3
// coefficient of restitution determines the elasticity of a collision: C_REST = [0,1]
//  if C_REST = 0 -> perfectly inelastic (particles stick together)
//  if C_REST = 1 -> perfectly elastic (no loss of speed)
#define C_REST 0.5
// set the iteration times
#define iteration_times 10000
// Must set 0 if run on Pi
#define NOT_RUN_ON_PI 1

struct body {
    double x, y; // position
    double vx, vy; // velocity
    double m; // mass
    double r; // radius of the particle
};


struct Quad {
	double x, y;
	double xlen, ylen;
};


void initQ(struct Quad *q, double a, double b, double c, double d){
	q->x = a;
	q->y = b;
	q->xlen = c;
	q->ylen = d;
};

bool contains(struct Quad *q, double a, double b) {
	if (a >= q->x && b >= q->y && a < q->x + q->xlen && b < q->y + q->ylen) return true;
	return false;
}


struct BHTree {
	struct body body;
	struct Quad quad;
	struct BHTree *NW;
	struct BHTree *NE;
	struct BHTree *SW;
	struct BHTree *SE;
};

void initB (struct BHTree *p, struct Quad q){
	p->body.x = p->body.y = 0;
	p->body.vx = p->body.vy = 0;
	p->body.m = p->body.r = 0;
	p->quad = q;
	p->NE = p->NW = p->SW = p->SE = NULL;
}

void insert_to_child(struct BHTree *p, struct body b);

void insert(struct BHTree *p, struct body b) {
	//printf("%f\n", b.x);
	struct body a = p->body;
	p->body.m += b.m;
	p->body.x = (p->body.x + b.x * b.m) / p->body.m;
	p->body.y = (p->body.y + b.y * b.m) / p->body.m;
	if (a.m == 0) {
	} else if (!(p->NW == NULL && p->NE == NULL && p->SW == NULL & p->SE == NULL)) {
		insert_to_child(p, b);
	} else {
		insert_to_child(p, a);
		insert_to_child(p, b);
	}
}

void insert_to_child(struct BHTree *p, struct body b){
	struct Quad *q = malloc(sizeof(struct Quad));
	if (b.x < p->quad.x + p->quad.xlen / 2) {
		if(b.y < p->quad.y + p->quad.ylen / 2) {
			if(p->SW == NULL) {
				p->SW = malloc(sizeof(struct BHTree));
				initQ(q, p->quad.x, p->quad.y, p->quad.xlen/2,
						p->quad.ylen/2); 
				initB(p->SW, *q);
			}
			insert(p->SW, b);
		} else {
			if(p->NW == NULL) {
				p->NW = malloc(sizeof(struct BHTree));
				initQ(q, p->quad.x, p->quad.y + p->quad.ylen/2, p->quad.xlen/2,
						p->quad.ylen/2); 
				initB(p->NW, *q);
			}
			insert(p->NW, b);
		}
	} else {
		if(b.y < p->quad.y + p->quad.ylen / 2) {
			if(p->SE == NULL) {
				p->SE = malloc(sizeof(struct BHTree));
				initQ(q, p->quad.x+p->quad.xlen/2, p->quad.y, p->quad.xlen/2,
						p->quad.ylen/2); 
				initB(p->SE, *q);
			}
			insert(p->SE, b);
		} else {
			if(p->NE == NULL) {
				p->NE = malloc(sizeof(struct BHTree));
				initQ(q, p->quad.x + p->quad.xlen/2, p->quad.y + p->quad.ylen/2, p->quad.xlen/2,
						p->quad.ylen/2); 
				initB(p->NE, *q);
			}
			/*printf("m is %f\n", b.m);
			printf("x, y is %f %f\n", b.x, b.y);
			printf("Select NE\n");
			*/
			insert(p->NE, b);
		}
	}
	free(q);
}

struct world {
    struct body *bodies;
	struct BHTree *tree;
    int num_bodies;
};

double max(double x, double y){if (x < y) return y; return x;}

void ooforce(struct body a, struct body b, double *fx, double *fy);

void computeForce(struct body b, struct BHTree *tree, double *fx, double *fy) {
	double s = max(tree->quad.xlen, tree->quad.ylen);
	double d = sqrt(pow(tree->body.x-b.x,2) + pow(tree->body.y-b.y,2));
	bool isleaf = tree->NW == NULL && tree->NE == NULL && tree->SW == NULL &&
		tree->SE == NULL;
	if (isleaf && tree->body.x == b.x && tree->body.y == b.y) {
		// same body
	} else if (s / d < THETA || isleaf){
		ooforce(b, tree->body, fx, fy);
	} else {
		if (tree->SE != NULL) computeForce(b, tree->SE, fx, fy);
		if (tree->NW != NULL) computeForce(b, tree->NW, fx, fy);
		if (tree->NE != NULL) computeForce(b, tree->NE, fx, fy);
		if (tree->SW != NULL) computeForce(b, tree->SW, fx, fy);
	}
}

void ooforce(struct body a, struct body b, double *fx, double *fy) {
	// Compute the x and y distances and total distance d between
	// bodies i and j
	double diff_x = b.x - a.x;
	double diff_y = b.y - a.y;
	double d = sqrt((diff_x * diff_x) + (diff_y * diff_y));

	if (d < 25) {
		d = 25;
	}
	double d_cubed = d * d * d;
	// Add force due to j to total force on i
	*(fx) += GRAV * (a.m * b.m / d_cubed) * diff_x;
	*(fy) += GRAV * (a.m * b.m / d_cubed) * diff_y;
}

clock_t total_time = 0;
//total_time.sec = 0;
//total_time.usec = 0;


/* This function initializes each particle's mass, velocity and position */
struct world* create_world(int num_bodies) {
    struct world *world = malloc(sizeof(struct world));

    world->num_bodies = num_bodies;
    world->bodies = malloc(sizeof(struct body)*num_bodies);

    int i = 0;
    double x;
    double y;
    double rc;

    int min_dim = (WIDTH < HEIGHT) ? WIDTH : HEIGHT;

    while (i < num_bodies / 100) {
        x = drand48() * WIDTH;
        y = drand48() * HEIGHT;
        rc = sqrt((WIDTH/2-x)*(WIDTH/2-x) + (y-HEIGHT/2)*(y-HEIGHT/2));
        if (rc <= min_dim/2) {
			int j;
			for (j = 0; j < 100; ++j){
				int k = i * 100 + j;
				world->bodies[k].x = x + drand48();
				world->bodies[k].y = y + drand48();

				world->bodies[k].vx = V_SCALAR * (y-HEIGHT/2) / rc;
				world->bodies[k].vy = V_SCALAR * (WIDTH/2-x) / rc;

				world->bodies[k].m = (1 / (0.025 + drand48())) * M_SCALAR;
				world->bodies[k].r = sqrt(world->bodies[k].m / M_PI) * R_SCALAR;
			}
            i++;
        }
    }

	
    return world;
}

// set the foreground color given RGB values between 0..255.
void set_color(Display *disp, GC gc, int r, int g, int b){
  unsigned long int p ;

  if (r < 0) r = 0; else if (r > 255) r = 255;
  if (g < 0) g = 0; else if (g > 255) g = 255;
  if (b < 0) b = 0; else if (b > 255) b = 255;

  p = (r << 16) | (g  << 8) | (b) ;

  XSetForeground(disp, gc, p) ;
}


/* This function updates the screen with the new positions of each particle */
void draw_world(Display *disp, Pixmap back_buf, GC gc, struct world *world) {
    int i;
    double x, y, r, r2;

    // we turn off aliasing for faster draws
    set_color(disp, gc, 255, 255, 255);
    XFillRectangle(disp, back_buf, gc, 0, 0, WIDTH, HEIGHT);

    for (i = 0; i < world->num_bodies; i++) {
        r = world->bodies[i].r;
        x = world->bodies[i].x - r;
        y = world->bodies[i].y - r;
        r2 = r + r;

        // draw body
        set_color(disp, gc, 255*7/10, 255*7/10, 255*7/10);
        XFillArc(disp, back_buf, gc, x, y, r2, r2, 0, 360*64); 
        set_color(disp, gc, 0, 0, 0);
        XDrawArc(disp, back_buf, gc, x, y, r2, r2, 0, 360*64); 
    }
}

void collision_step(struct world *world) {
    int a, b;
    double r, x, y, vx, vy;

    // Impose screen boundaries by reversing direction if body is off screen
    for (a = 0; a < world->num_bodies; a++) {
        r = world->bodies[a].r;
        x = world->bodies[a].x;
        y = world->bodies[a].y;
        vx = world->bodies[a].vx;
        vy = world->bodies[a].vy;

        if (x-r < 0) { // left edge
            if (vx < 0) { world->bodies[a].vx = -C_REST * vx; }
            world->bodies[a].x = r;
        } else if (x+r > WIDTH) { // right edge
            if (vx > 0) { world->bodies[a].vx = -C_REST * vx; }
            world->bodies[a].x = WIDTH - r;
        }

        if (y-r < 0) { // bottom edge
            if (vy < 0) { world->bodies[a].vy = -C_REST * vy; }
            world->bodies[a].y = r;
        } else if (y+r > HEIGHT) { // top edge
            if (vy > 0) { world->bodies[a].vy = -C_REST * vy; }
            world->bodies[a].y = HEIGHT - r;
        }
    }
}

void position_step(struct world *world, double time_res, int thread_count) {
    int i, j;
    double d, d_cubed, diff_x, diff_y;


    /* The forces array stores the x and y components of the total force acting
     * on each body. The forces are index like this:
     *     F on body i in the x dir = F_x[i]
     *     F on body i in the y dir = F_y[i] */
    double *force_x = (double*)malloc(sizeof(double) * world->num_bodies);
	double *force_y = (double*)malloc(sizeof(double) * world->num_bodies);
    // initialize all forces to zero
    force_x = memset(force_x, 0, sizeof(double) * world->num_bodies);
	force_y = memset(force_y, 0, sizeof(double) * world->num_bodies);

	struct Quad *q = malloc(sizeof(struct Quad));
	initQ(q, 0, 0, WIDTH, HEIGHT);
	world->tree = malloc(sizeof(struct BHTree));
	initB(world->tree, *q);
	for (i = 0; i < world->num_bodies; ++i){
		insert(world->tree, world->bodies[i]);
	}
    /* Compute the net force on each body */
    /*
	for (i = 0; i < world->num_bodies; i++) {
        for (j = 0; j < world->num_bodies; j++) {
            if (i == j) {
                continue;
            }
            // Compute the x and y distances and total distance d between
            // bodies i and j
            diff_x = world->bodies[j].x - world->bodies[i].x;
            diff_y = world->bodies[j].y - world->bodies[i].y;
            d = sqrt((diff_x * diff_x) + (diff_y * diff_y));

			if (d < 25) {
                d = 25;
            }
            d_cubed = d * d * d;
            // Add force due to j to total force on i
            force_x[i] += GRAV * (world->bodies[i].m * world->bodies[j].m
                    / d_cubed) * diff_x;
            force_y[i] += GRAV * (world->bodies[i].m * world->bodies[j].m
                    / d_cubed) * diff_y;
        }
    }*/

	# pragma omp parallel num_threads(thread_count) 
	{
		int my_rank = omp_get_thread_num();
		int thread_count = omp_get_num_threads();

		int local_bodies = world->num_bodies/thread_count;
		int local_index = local_bodies*my_rank;
		int elem;

		for (i = 0; i < local_bodies; ++i) {
			computeForce(world->bodies[i+local_index], world->tree,
					&force_x[i+local_index], &force_y[i+local_index]);
		}

		// Update the velocity and position of each body
		
		for (i = 0; i < local_bodies; i++) {
			elem = i + local_index;
			// Update velocities
			world->bodies[elem].vx += force_x[elem] * time_res / world->bodies[elem].m;
			world->bodies[elem].vy += force_y[elem] * time_res / world->bodies[elem].m;		
			
			// Update positions
			world->bodies[elem].x += world->bodies[elem].vx * time_res;
			world->bodies[elem].y += world->bodies[elem].vy * time_res;
		}
	}

}

void step_world(struct world *world, double time_res, int thread_count) {
	
	struct tms ttt;
	clock_t start, end;
	start = times(&ttt);
    position_step(world, time_res, thread_count);
	end = times(&ttt);
	total_time += end - start;

    collision_step(world);
}


/* Main method runs initialize() and update() */
int main(int argc, char **argv) {
	//total_time.tv_sec = 0;
	//total_time.tv_usec = 0;
    /* get num bodies from the command line */
    int num_bodies;
    num_bodies = (argc == 3) ? atoi(argv[1]) : DEF_NUM_BODIES;
    int thread_count = (argc == 3) ? atoi(argv[2]):1;
	printf("Universe has %d bodies.\n", num_bodies);
	printf("%d threads \n", thread_count);

    /* set up the universe */
    time_t cur_time;
    time(&cur_time);
    srand48((long)cur_time); // seed the RNG used in create_world
    struct world *world = create_world(num_bodies);
    /* set up graphics using Xlib */
#if NOT_RUN_ON_PI
    Display *disp = XOpenDisplay(NULL);
    int scr = DefaultScreen(disp);
    Window win = XCreateSimpleWindow(
            disp,
            RootWindow(disp, scr),
            0, 0,
            WIDTH, HEIGHT,
            0,
            BlackPixel(disp, scr), WhitePixel(disp, scr));
    XStoreName(disp, win, "N-Body Simulator");

    Pixmap back_buf = XCreatePixmap(disp, RootWindow(disp, scr),
            WIDTH, HEIGHT, DefaultDepth(disp, scr));
    GC gc = XCreateGC(disp, back_buf, 0, 0);

    // Make sure we're only looking for messages about closing the window
    Atom del_window = XInternAtom(disp, "WM_DELETE_WINDOW", 0);
    XSetWMProtocols(disp, win, &del_window, 1);

    XSelectInput(disp, win, StructureNotifyMask);
    XMapWindow(disp, win);
    XEvent event;
    // wait until window is mapped
    while (1) {
        XNextEvent(disp, &event);
        if (event.type == MapNotify) {
            break;
        }
    }
#endif

    struct timespec delay={0, 1000000000 / 60}; // for 60 FPS
    struct timespec remaining;
	double delta_t = 0.1;
	int ii;
	
    for(ii = 0; ii < iteration_times; ii++){
        // check if the window has been closed
#if NOT_RUN_ON_PI	
        if (XCheckTypedEvent(disp, ClientMessage, &event)) {
            break;
        }

        // we first draw to the back buffer then copy it to the front (`win`)
        draw_world(disp, back_buf, gc, world);
        XCopyArea(disp, back_buf, win, gc, 0, 0, WIDTH, HEIGHT, 0, 0);
#endif

        step_world(world, delta_t, thread_count);

		//if you want to watch the process in 60 FPS
        //nanosleep(&delay, &remaining);
    }

//	printf("Total Time = %f\n", (double)total_time.tv_sec + (double)total_time.tv_usec/1000000);
	printf("Nbody Position Calculation Time = :%lf s\n",(double)total_time / (sysconf(_SC_CLK_TCK)));

#if NOT_RUN_ON_PI	
    XFreeGC(disp, gc);
    XFreePixmap(disp, back_buf);
    XDestroyWindow(disp, win);
    XCloseDisplay(disp);
#endif

    return 0;
}
