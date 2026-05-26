document.addEventListener("DOMContentLoaded", () => {

  const canvas = document.getElementById("bg-canvas");

  if (!canvas) return;

  const ctx = canvas.getContext("2d");

  let width;
  let height;
  let particles = [];

  const COLORS = [
    [168, 85, 247],
    [124, 58, 237],
    [196, 148, 255],
    [147, 51, 234]
  ];

  function resize() {
    width = canvas.width = window.innerWidth;
    height = canvas.height = window.innerHeight;
  }

  function random(min, max) {
    return min + Math.random() * (max - min);
  }

  function createParticle() {

    const color =
      COLORS[Math.floor(Math.random() * COLORS.length)];

    return {
      x: random(0, width),
      y: random(0, height),

      vx: random(-0.15, 0.15),
      vy: random(-0.1, 0.1),

      radius: random(1.2, 2.6),
      alpha: random(0.25, 0.65),

      color
    };
  }

  function init() {

    resize();

    particles = [];

    const count =
      Math.min(Math.floor((width * height) / 16000), 80);

    for (let i = 0; i < count; i++) {
      particles.push(createParticle());
    }
  }

  const mouse = {
    x: -999,
    y: -999
  };

  window.addEventListener("mousemove", (e) => {
    mouse.x = e.clientX;
    mouse.y = e.clientY;
  });

  function animate() {

    ctx.clearRect(0, 0, width, height);

    const LINK_DISTANCE = 130;

    for (let i = 0; i < particles.length; i++) {

      const a = particles[i];

      for (let j = i + 1; j < particles.length; j++) {

        const b = particles[j];

        const dx = a.x - b.x;
        const dy = a.y - b.y;

        const distance =
          Math.sqrt(dx * dx + dy * dy);

        if (distance < LINK_DISTANCE) {

          ctx.beginPath();

          ctx.moveTo(a.x, a.y);
          ctx.lineTo(b.x, b.y);

          const opacity =
            (1 - distance / LINK_DISTANCE) * 0.18;

          ctx.strokeStyle =
            `rgba(168,85,247,${opacity})`;

          ctx.lineWidth = 0.8;

          ctx.stroke();
        }
      }

      const mdx = a.x - mouse.x;
      const mdy = a.y - mouse.y;

      const mouseDistance =
        Math.sqrt(mdx * mdx + mdy * mdy);

      if (mouseDistance < 180) {

        ctx.beginPath();

        ctx.moveTo(a.x, a.y);
        ctx.lineTo(mouse.x, mouse.y);

        const opacity =
          (1 - mouseDistance / 180) * 0.18;

        ctx.strokeStyle =
          `rgba(${a.color.join(",")},${opacity})`;

        ctx.lineWidth = 0.8;

        ctx.stroke();
      }

      ctx.beginPath();

      ctx.arc(
        a.x,
        a.y,
        a.radius,
        0,
        Math.PI * 2
      );

      ctx.fillStyle =
        `rgba(${a.color.join(",")},${a.alpha})`;

      ctx.fill();

      a.x += a.vx;
      a.y += a.vy;

      if (a.x < 0) a.x = width;
      if (a.x > width) a.x = 0;

      if (a.y < 0) a.y = height;
      if (a.y > height) a.y = 0;
    }

    requestAnimationFrame(animate);
  }

  window.addEventListener("resize", init);

  init();

  animate();

});