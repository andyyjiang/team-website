// Footer year
document.getElementById("year").textContent = new Date().getFullYear();

// Mobile menu toggle
const toggle = document.querySelector(".nav-toggle");
const links = document.getElementById("nav-links");

toggle.addEventListener("click", () => {
  const open = links.classList.toggle("open");
  toggle.setAttribute("aria-expanded", String(open));
});

// Close the menu after picking a link
links.addEventListener("click", (event) => {
  if (event.target.tagName === "A") {
    links.classList.remove("open");
    toggle.setAttribute("aria-expanded", "false");
  }
});
