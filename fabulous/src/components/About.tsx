import React from 'react';
import { useState } from 'react';
import "./About.css";
import Stephen from "../assets/stephen.png";
import Linh from "../assets/linh.png";
import Benjamin from "../assets/benjamin.png";
import Medy from "../assets/medy.png";


interface TeamMember {
  name: string;
  photo: string;
  bio: string;
}

const teamMembers: TeamMember[] = [
  {
    name: 'Stephen',
    photo: Stephen,
    bio: 'Stephen is our multiculti project manager and also part of the design team. He is a great team player and always brings a positive attitude to the table. Also making sure that the project is on track and that the team is happy.',
  },
  {
    name: 'Benjamin',
    photo: Benjamin,
    bio: 'Benjamin was at the root of the project. With his scientific background, he now oversees the technical aspects of the project. He is also part of the developping team and is always looking for ways to improve the project.',
  },
  {
    name: 'Linh',
    photo: Linh,
    bio: "Linh has an outstanding talent to make things feel beautiful. At the heart of our design team, she's always looking for ways to make the project more user-friendly and visually appealing.",
  },
  {
    name: 'Medy',
    photo: Medy,
    bio: "A strong developper lies in his core. Medy is the person you want to speak with when it comes to making things work. Beside integrating all the designer's dreams and the scientists wishes, he also is responsible for everything you see working flawlessly.",
  },
];

const About: React.FC = () => {
  const [isOpen, setIsOpen] = useState(false);

  const handleOpen = () => {
    setIsOpen(true);
  };

  const handleClose = () => {
    setIsOpen(false);
  };

  return (
    <div>
      <a onClick={handleOpen}>About</a>
      {isOpen && (
        <div className="overlay">
          <div className="overlay-content">
            <button className="close" onClick={handleClose}>
              ×
            </button>
            <h2>About Our Team</h2>
            <ul>
              {teamMembers.map((member) => (
                <li key={member.name}>
                  <img src={member.photo} alt={member.name} className='pics'/>
                  <h3>{member.name}</h3>
                  <p>{member.bio}</p>
                </li>
              ))}
            </ul>
          </div>
        </div>
      )}
    </div>
  );
};

export default About;