document.addEventListener('DOMContentLoaded', function() {

  // Collapsible sections
  var coll = document.getElementsByClassName("collapsible");
  for (var i = 0; i < coll.length; i++) {
    coll[i].addEventListener("click", function() {
      this.classList.toggle("active");
      var content = this.nextElementSibling;
      if (content.style.maxHeight) {
        content.style.maxHeight = null;
      } else {
        content.style.maxHeight = (content.scrollHeight + 2) + "px";
      }
    });
  }

  // Back to top button
  var mybutton = document.getElementById("topBtn");
  window.onscroll = function() {
    if (document.body.scrollTop > 20 || document.documentElement.scrollTop > 20) {
      mybutton.style.display = "block";
    } else {
      mybutton.style.display = "none";
    }
  };

});

// These stay outside as they're called directly from onclick attributes
function topFunction() {
  window.scrollTo({ top: 0, behavior: 'smooth' });
}

function myFunction() {
  var x = document.getElementById("myTopnav");
  if (x.className === "topnav") {
    x.className += " responsive";
  } else {
    x.className = "topnav";
  }
}

function scrollToCard(id) {
  document.querySelectorAll('.card').forEach(c => c.classList.remove('highlighted'));
  const el = document.getElementById(id);
  el.scrollIntoView({ behavior: 'smooth', block: 'center' });
  el.classList.add('highlighted');
  setTimeout(() => el.classList.remove('highlighted'), 1500);
}